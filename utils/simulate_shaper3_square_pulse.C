#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int kEvents = 10000;
constexpr int kPulsesPerEvent = 400;
constexpr double kCurrentUa = 0.4;
constexpr double kPulseDurationNs = 1.0;
constexpr double kMaxDelayNs = 100.0;
constexpr double kPositionPerNsUm = 50.0;
constexpr double kThresholdMv = 1.0;
constexpr double kSampleStepNs = 0.1;

std::string filenameNumber(double value) {
  std::ostringstream formatted;
  formatted << std::fixed << std::setprecision(2) << value;
  std::string result = formatted.str();
  while (!result.empty() && result.back() == '0') {
    result.pop_back();
  }
  if (!result.empty() && result.back() == '.') {
    result.pop_back();
  }
  std::replace(result.begin(), result.end(), '.', 'p');
  return result;
}

struct ShaperState {
  double stage1 = 0.0;
  double stage2 = 0.0;
  double stage3 = 0.0;
  double stage4 = 0.0;
};

ShaperState derivative(const ShaperState& state, double inputFcPerNs,
                       double polePerNs) {
  return {-polePerNs * state.stage1 + inputFcPerNs,
          state.stage1 - polePerNs * state.stage2,
          state.stage2 - polePerNs * state.stage3,
          state.stage3 - polePerNs * state.stage4};
}

ShaperState addScaled(const ShaperState& state, const ShaperState& increment,
                      double scale) {
  return {state.stage1 + scale * increment.stage1,
          state.stage2 + scale * increment.stage2,
          state.stage3 + scale * increment.stage3,
          state.stage4 + scale * increment.stage4};
}

void advanceShaper(ShaperState& state, double inputFcPerNs,
                   double shapingTimeNs, double stepNs) {
  const double polePerNs = 3.0 / shapingTimeNs;
  const ShaperState k1 = derivative(state, inputFcPerNs, polePerNs);
  const ShaperState k2 =
      derivative(addScaled(state, k1, stepNs / 2.0), inputFcPerNs, polePerNs);
  const ShaperState k3 =
      derivative(addScaled(state, k2, stepNs / 2.0), inputFcPerNs, polePerNs);
  const ShaperState k4 =
      derivative(addScaled(state, k3, stepNs), inputFcPerNs, polePerNs);
  state.stage1 += stepNs / 6.0 *
                  (k1.stage1 + 2.0 * k2.stage1 + 2.0 * k3.stage1 + k4.stage1);
  state.stage2 += stepNs / 6.0 *
                  (k1.stage2 + 2.0 * k2.stage2 + 2.0 * k3.stage2 + k4.stage2);
  state.stage3 += stepNs / 6.0 *
                  (k1.stage3 + 2.0 * k2.stage3 + 2.0 * k3.stage3 + k4.stage3);
  state.stage4 += stepNs / 6.0 *
                  (k1.stage4 + 2.0 * k2.stage4 + 2.0 * k3.stage4 + k4.stage4);
}

double shaperOutputMv(const ShaperState& state,
                      double gainScaledNormalization) {
  return gainScaledNormalization * state.stage4;
}

}  // namespace

void simulate_shaper3_square_pulse(
    double shapingTimeNs = 25.0, double gainMvPerFc = 2.0,
    double timeSpreadNs = 100.0) {
  if (!std::isfinite(shapingTimeNs) || shapingTimeNs <= 0.0) {
    std::cerr << "The shaping time must be a positive finite value."
              << std::endl;
    return;
  }
  if (!std::isfinite(gainMvPerFc) || gainMvPerFc <= 0.0) {
    std::cerr << "The gain must be a positive finite value." << std::endl;
    return;
  }
  if (!std::isfinite(timeSpreadNs) || timeSpreadNs < 0.0) {
    std::cerr << "TimeSpread must be a non-negative finite value." << std::endl;
    return;
  }
  const std::string outputBase =
      "plot/DTrelation_tau" + filenameNumber(shapingTimeNs) + "_spread" +
      filenameNumber(timeSpreadNs);
  const std::string plotPath = outputBase + ".png";
  const std::string rootPath = outputBase + ".root";

  TGraph riseTimeGraph;
  riseTimeGraph.SetName("rise_time_vs_position");
  riseTimeGraph.SetTitle(
      "Third-order shaper: 1 mV threshold crossing vs position;"
      "Position [#mum];Threshold crossing time [ns]");
  riseTimeGraph.SetMarkerStyle(20);
  riseTimeGraph.SetMarkerColor(kBlue + 1);
  riseTimeGraph.SetMarkerSize(0.45);

  auto outputFile =
      std::unique_ptr<TFile>(TFile::Open(rootPath.c_str(), "RECREATE"));
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Cannot create ROOT output file: " << rootPath << std::endl;
    return;
  }

  TTree events("events", "Square-pulse shaper simulation events");
  ULong64_t event = 0;
  double delayNs = 0.0;
  double positionUm = 0.0;
  double riseTimeNs = -1.0;
  bool riseThresholdCrossed = false;
  events.Branch("event", &event, "event/l");
  events.Branch("delay_ns", &delayNs, "delay_ns/D");
  events.Branch("position_um", &positionUm, "position_um/D");
  events.Branch("t_rise_ns", &riseTimeNs, "t_rise_ns/D");
  events.Branch("rise_threshold_crossed", &riseThresholdCrossed,
                "rise_threshold_crossed/O");

  TRandom3 random(12345);
  const double gainScaledNormalization =
      gainMvPerFc * 6.0 * std::exp(3.0) / std::pow(shapingTimeNs, 3);
  TGraph waveform;
  waveform.SetName("shaper3_waveform_event0");
  waveform.SetTitle(
      "Third-order shaper response to a 400-pulse train;Time [ns];"
      "Output voltage [mV]");
  int eventsWithoutCrossing = 0;
  for (event = 0; event < kEvents; ++event) {
    delayNs = random.Uniform(0.0, kMaxDelayNs);
    positionUm = delayNs * kPositionPerNsUm;

    const double latestPulseEndNs = delayNs + timeSpreadNs + kPulseDurationNs;
    const double waveformEndNs = latestPulseEndNs + 8.0 * shapingTimeNs;
    const int waveformSamples =
        static_cast<int>(std::ceil(waveformEndNs / kSampleStepNs)) + 1;
    std::vector<double> inputFcPerNs(waveformSamples, 0.0);
    for (int pulse = 0; pulse < kPulsesPerEvent; ++pulse) {
      const double pulseStartNs =
          random.Uniform(delayNs, delayNs + timeSpreadNs);
      const double pulseEndNs = pulseStartNs + kPulseDurationNs;
      const int firstBin =
          static_cast<int>(std::floor(pulseStartNs / kSampleStepNs));
      const int endBin =
          static_cast<int>(std::ceil(pulseEndNs / kSampleStepNs));
      for (int bin = firstBin; bin < endBin && bin < waveformSamples; ++bin) {
        const double binStartNs = bin * kSampleStepNs;
        const double overlapNs =
            std::max(0.0, std::min(pulseEndNs, binStartNs + kSampleStepNs) -
                              std::max(pulseStartNs, binStartNs));
        inputFcPerNs[bin] += kCurrentUa * overlapNs / kSampleStepNs;
      }
    }

    ShaperState state;
    double previousOutputMv = 0.0;
    double previousTimeNs = 0.0;
    riseTimeNs = -1.0;
    riseThresholdCrossed = false;
    for (int sample = 0; sample < waveformSamples; ++sample) {
      const double timeNs = sample * kSampleStepNs;
      const double outputMv = shaperOutputMv(state, gainScaledNormalization);
      if (!riseThresholdCrossed && sample > 0 &&
          previousOutputMv < kThresholdMv &&
          outputMv >= kThresholdMv) {
        const double fraction =
            (kThresholdMv - previousOutputMv) /
            (outputMv - previousOutputMv);
        riseTimeNs = previousTimeNs +
                     fraction * (timeNs - previousTimeNs);
        riseThresholdCrossed = true;
      }
      if (event == 0) {
        waveform.SetPoint(sample, timeNs, outputMv);
      }
      if (!riseThresholdCrossed) {
        previousOutputMv = outputMv;
        previousTimeNs = timeNs;
      }
      if (sample + 1 < waveformSamples) {
        advanceShaper(state, inputFcPerNs[sample], shapingTimeNs,
                      kSampleStepNs);
      }
    }

    if (riseThresholdCrossed) {
      riseTimeGraph.SetPoint(riseTimeGraph.GetN(), positionUm, riseTimeNs);
    } else {
      ++eventsWithoutCrossing;
    }
    events.Fill();
  }

  TCanvas canvas("c_shaper3_rise_time",
                 "Third-order shaper rise time versus position", 1000, 700);
  riseTimeGraph.Draw("AP");
  canvas.SaveAs(plotPath.c_str());

  outputFile->cd();
  riseTimeGraph.Write();
  waveform.Write();
  events.Write();
  outputFile->Close();

  std::ifstream plotFile(plotPath, std::ios::binary);
  if (!plotFile) {
    std::cerr << "Failed to create rise-time plot: " << plotPath << std::endl;
    return;
  }
  std::cout << "Simulated " << kEvents << " events with " << kPulsesPerEvent
            << " pulses per event and TimeSpread = " << timeSpreadNs
            << " ns. Events without a 1 mV crossing: " << eventsWithoutCrossing
            << ". Plot saved to " << plotPath << "; event data saved to "
            << rootPath << std::endl;
}
