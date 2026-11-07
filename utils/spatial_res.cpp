#include <TAxis.h>
#include <array>
#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TSpline.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

namespace {

constexpr int kSlices = 23;
constexpr double kSplineMinimumTrackY = -0.65;
constexpr double kSplineMaximumTrackY = 0.0;
constexpr int kResidualBins = 7;

double crystalBall(double* x, double* parameters) {
  const double alpha = std::abs(parameters[3]);
  const double exponent = parameters[4];
  if (alpha <= 0.0 || exponent <= 0.0 || parameters[2] <= 0.0) {
    return 0.0;
  }

  const double standardized = (x[0] - parameters[1]) / parameters[2];
  if (standardized > -alpha) {
    return parameters[0] * std::exp(-0.5 * standardized * standardized);
  }

  const double coefficient =
      std::pow(exponent / alpha, exponent) * std::exp(-0.5 * alpha * alpha);
  const double offset = exponent / alpha - alpha;
  return parameters[0] * coefficient *
         std::pow(offset - standardized, -exponent);
}

std::string twoDigits(int value) {
  std::ostringstream result;
  result << std::setw(2) << std::setfill('0') << value;
  return result.str();
}

bool invertSpline(const TSpline3& spline, double minimum, double maximum,
                  double time, double& trackY) {
  constexpr int kSearchSteps = 2048;
  constexpr double kTimeTolerance = 1.e-10;
  constexpr double kPositionTolerance = 1.e-8;
  int solutions = 0;
  double solution = 0.0;
  double left = minimum;
  double leftValue = spline.Eval(left) - time;

  const auto addSolution = [&](double candidate) {
    if (solutions == 0 ||
        std::abs(candidate - solution) > kPositionTolerance) {
      solution = candidate;
      ++solutions;
    }
  };

  if (std::abs(leftValue) <= kTimeTolerance) {
    addSolution(left);
  }
  for (int step = 1; step <= kSearchSteps; ++step) {
    const double right =
        minimum + (maximum - minimum) * step / kSearchSteps;
    const double rightValue = spline.Eval(right) - time;
    if (std::abs(rightValue) <= kTimeTolerance) {
      addSolution(right);
    } else if ((leftValue < 0.0 && rightValue > 0.0) ||
               (leftValue > 0.0 && rightValue < 0.0)) {
      double low = left;
      double high = right;
      double lowValue = leftValue;
      for (int iteration = 0; iteration < 60; ++iteration) {
        const double middle = 0.5 * (low + high);
        const double middleValue = spline.Eval(middle) - time;
        if (std::abs(middleValue) <= kTimeTolerance) {
          low = middle;
          high = middle;
          break;
        }
        if ((lowValue < 0.0 && middleValue > 0.0) ||
            (lowValue > 0.0 && middleValue < 0.0)) {
          high = middle;
        } else {
          low = middle;
          lowValue = middleValue;
        }
      }
      addSolution(0.5 * (low + high));
    }
    left = right;
    leftValue = rightValue;
  }

  if (solutions != 1) {
    return false;
  }
  trackY = solution;
  return true;
}

}  // namespace

void spatial_res(const char* runNumber) {
  if (!runNumber || runNumber[0] == '\0' ||
      !std::all_of(runNumber, runNumber + std::char_traits<char>::length(runNumber),
                   [](unsigned char character) {
                     return character >= '0' && character <= '9';
                   })) {
    std::cerr << "A numeric run number is required." << std::endl;
    return;
  }

  const std::string runBase = std::string("results/run") + runNumber;
  const std::string inputPath = runBase + "_correlations.root";
  const std::string outputDirectory = runBase + "_spatial_res";
  const std::string plotPath = runBase + "_spatial_res.png";
  const std::string rootPath = runBase + "_spatial_res.root";

  std::unique_ptr<TFile> inputFile(TFile::Open(inputPath.c_str(), "READ"));
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Cannot open input file: " << inputPath << std::endl;
    return;
  }

  auto* correlation =
      dynamic_cast<TH2*>(inputFile->Get("h2_riseTime_vs_trackY"));
  if (!correlation) {
    correlation = dynamic_cast<TH2*>(inputFile->Get("riseTime_vs_trackY"));
  }
  if (!correlation) {
    std::cerr << "Histogram 'h2_riseTime_vs_trackY' not found in "
              << inputPath << std::endl;
    return;
  }
  if (correlation->GetNbinsX() < kSlices) {
    std::cerr << "The input histogram has fewer than " << kSlices
              << " bins on the trackY axis." << std::endl;
    return;
  }

  if (gSystem->mkdir(outputDirectory.c_str(), true) != 0 &&
      gSystem->AccessPathName(outputDirectory.c_str())) {
    std::cerr << "Cannot create output directory: " << outputDirectory
              << std::endl;
    return;
  }

  const int sourceBins = correlation->GetNbinsX();
  const double xMinimum = correlation->GetXaxis()->GetXmin();
  const double xMaximum = correlation->GetXaxis()->GetXmax();
  TGraphErrors mpvGraph;
  mpvGraph.SetName("crystal_ball_mpv_vs_trackY");
  mpvGraph.SetTitle(
      "Crystal Ball MPV vs track y;Track y [cm];Time-at-threshold [ns]");
  mpvGraph.SetMarkerStyle(20);
  mpvGraph.SetMarkerColor(kRed + 1);
  mpvGraph.SetLineColor(kRed + 1);
  mpvGraph.SetMarkerSize(1.0);
  TGraphErrors timeResolutionGraph;
  timeResolutionGraph.SetName("crystal_ball_sigma_vs_trackY");
  timeResolutionGraph.SetTitle(
      "Time resolution vs track y;Track y slice center [cm];"
      "Crystal Ball #sigma [ns]");
  timeResolutionGraph.SetMarkerStyle(20);
  timeResolutionGraph.SetMarkerColor(kBlue + 1);
  timeResolutionGraph.SetLineColor(kBlue + 1);
  timeResolutionGraph.SetMarkerSize(1.0);
  TGraph timeRmsGraph;
  timeRmsGraph.SetName("time_at_threshold_rms_vs_trackY");
  timeRmsGraph.SetTitle(
      "Time distribution RMS vs track y;Track y slice center [cm];"
      "Time width [ns]");
  timeRmsGraph.SetMarkerStyle(21);
  timeRmsGraph.SetMarkerColor(kRed + 1);
  timeRmsGraph.SetLineColor(kRed + 1);

  for (int slice = 0; slice < kSlices; ++slice) {
    const int firstXBin = 1 + slice * sourceBins / kSlices;
    const int lastXBin = (slice + 1) * sourceBins / kSlices;
    const double xLow =
        xMinimum + (xMaximum - xMinimum) * slice / kSlices;
    const double xHigh =
        xMinimum + (xMaximum - xMinimum) * (slice + 1) / kSlices;
    const double xCenter = 0.5 * (xLow + xHigh);
    const std::string sliceId = twoDigits(slice + 1);
    const std::string histogramName =
        "time_at_threshold_slice_" + sliceId;
    const std::string sliceBase =
        outputDirectory + "/run" + runNumber + "_slice_" + sliceId;
    const std::string sliceRootPath = sliceBase + ".root";
    const std::string slicePlotPath = sliceBase + ".png";

    std::unique_ptr<TH1D> distribution(
        correlation->ProjectionY(histogramName.c_str(), firstXBin, lastXBin));
    if (!distribution) {
      std::cerr << "Failed to project trackY slice " << sliceId << "."
                << std::endl;
      return;
    }
    distribution->SetDirectory(nullptr);
    distribution->SetTitle(
        ("Time-at-threshold distribution, track y [" +
         std::to_string(xLow) + ", " + std::to_string(xHigh) +
         "];Time-at-threshold [ns];Events")
            .c_str());
    if (distribution->GetEntries() > 0.0) {
      timeRmsGraph.SetPoint(timeRmsGraph.GetN(), xCenter,
                            distribution->GetRMS());
    }

    std::unique_ptr<TFile> sliceFile(
        TFile::Open(sliceRootPath.c_str(), "RECREATE"));
    if (!sliceFile || sliceFile->IsZombie()) {
      std::cerr << "Cannot create slice file: " << sliceRootPath << std::endl;
      return;
    }

    TCanvas sliceCanvas(("c_slice_" + sliceId).c_str(),
                        "Time-at-threshold slice", 900, 650);
    distribution->SetLineColor(kBlue + 1);
    distribution->SetMarkerStyle(20);
    distribution->Draw("E");

    TLegend fitLegend(0.62, 0.78, 0.88, 0.88);
    bool fitSucceeded = false;
    if (distribution->GetEntries() > 0.0 && distribution->GetMaximum() > 0.0) {
      const double binWidth = distribution->GetXaxis()->GetBinWidth(1);
      const double sigmaStart =
          std::max(distribution->GetRMS(), binWidth);
      const double meanStart =
          distribution->GetXaxis()->GetBinCenter(distribution->GetMaximumBin());
      TF1 crystalBallFit(("crystal_ball_fit_" + sliceId).c_str(),
                         crystalBall, distribution->GetXaxis()->GetXmin(),
                         distribution->GetXaxis()->GetXmax(), 5);
      crystalBallFit.SetParNames("Amplitude", "MPV", "Sigma", "Alpha", "N");
      crystalBallFit.SetParameters(distribution->GetMaximum(), meanStart,
                                   sigmaStart, 1.5, 3.0);
      crystalBallFit.SetParLimits(0, 0.0, distribution->GetMaximum() * 100.0);
      crystalBallFit.SetParLimits(1, distribution->GetXaxis()->GetXmin(),
                                  distribution->GetXaxis()->GetXmax());
      crystalBallFit.SetParLimits(2, binWidth * 0.1,
                                  distribution->GetXaxis()->GetXmax() -
                                      distribution->GetXaxis()->GetXmin());
      crystalBallFit.SetParLimits(3, 0.1, 10.0);
      crystalBallFit.SetParLimits(4, 1.01, 100.0);
      const int fitStatus = distribution->Fit(&crystalBallFit, "QRS0");
      fitSucceeded = fitStatus == 0;
      if (fitSucceeded) {
        auto* fittedFunction =
            distribution->GetFunction(crystalBallFit.GetName());
        if (!fittedFunction) {
          std::cerr << "Fitted Crystal Ball function is missing for slice "
                    << sliceId << "." << std::endl;
          fitSucceeded = false;
        } else {
          fittedFunction->SetLineColor(kRed + 1);
          fittedFunction->SetLineWidth(3);
          fittedFunction->Draw("SAME");
          fitLegend.AddEntry(fittedFunction, "Crystal Ball fit", "l");
        }
        fitLegend.Draw();
        if (fitSucceeded) {
          mpvGraph.AddPoint(xCenter, crystalBallFit.GetParameter(1));
          mpvGraph.SetPointError(mpvGraph.GetN() - 1, 0.0,
                                 crystalBallFit.GetParError(1));
          timeResolutionGraph.AddPoint(xCenter,
                                       crystalBallFit.GetParameter(2));
          timeResolutionGraph.SetPointError(
              timeResolutionGraph.GetN() - 1, 0.0,
              crystalBallFit.GetParError(2));
        }
      } else {
        std::cerr << "Crystal Ball fit failed for slice " << sliceId
                  << " (status " << fitStatus << ")." << std::endl;
      }
    } else {
      std::cerr << "Not enough entries to fit slice " << sliceId << "."
                << std::endl;
    }

    sliceFile->cd();
    if (distribution->Write() <= 0 || sliceCanvas.Write() <= 0) {
      std::cerr << "Failed to write slice " << sliceId << " to "
                << sliceRootPath << std::endl;
      return;
    }
    sliceCanvas.SaveAs(slicePlotPath.c_str());
    sliceFile->Close();
    std::ifstream slicePlotFile(slicePlotPath, std::ios::binary);
    if (!slicePlotFile) {
      std::cerr << "Failed to create slice plot: " << slicePlotPath
                << std::endl;
      return;
    }
    std::cout << "Saved slice " << sliceId << " distribution to "
              << sliceRootPath << " and " << slicePlotPath
              << (fitSucceeded ? " (fit successful)." : " (fit unavailable).")
              << std::endl;
  }

  std::unique_ptr<TFile> outputFile(TFile::Open(rootPath.c_str(), "RECREATE"));
  if (!outputFile || outputFile->IsZombie()) {
    std::cerr << "Cannot create output file: " << rootPath << std::endl;
    return;
  }

  outputFile->cd();
  TCanvas timeResolutionCanvas("c_time_resolution",
                               "Time resolution vs track y", 1000, 700);
  TMultiGraph timeResolution;
  timeResolution.SetName("time_resolution_and_rms_vs_trackY");
  timeResolution.SetTitle(
      "Time resolution vs track y;Track y slice center [cm];Time width [ns]");
  timeResolution.Add(&timeResolutionGraph, "P");
  timeResolution.Add(&timeRmsGraph, "P");
  timeResolution.Draw("A");
  TLegend timeResolutionLegend(0.65, 0.76, 0.88, 0.89);
  timeResolutionLegend.AddEntry(&timeResolutionGraph, "Crystal Ball #sigma",
                                "p");
  timeResolutionLegend.AddEntry(&timeRmsGraph, "Distribution RMS", "p");
  timeResolutionLegend.Draw();
  if (timeResolutionGraph.GetN() == 0) {
    std::cerr << "No successful Crystal Ball fits; time-resolution graph "
                 "contains no points."
              << std::endl;
  }
  if (timeResolutionGraph.Write() <= 0 || timeRmsGraph.Write() <= 0 ||
      timeResolution.Write() <= 0 ||
      timeResolutionCanvas.Write() <= 0) {
    std::cerr << "Failed to write time-resolution output to " << rootPath
              << std::endl;
    return;
  }
  const std::string timeResolutionPlotPath = runBase + "_time_res.png";
  timeResolutionCanvas.SaveAs(timeResolutionPlotPath.c_str());

  TCanvas canvas("c_spatial_res", "Time-at-threshold vs track y", 1100, 800);
  correlation->SetTitle(
      "Time-at-threshold vs track y;Track y [cm];Time-at-threshold [ns]");
  correlation->Draw("COLZ");
  std::unique_ptr<TSpline3> mpvSpline;
  TGraph splinePoints;
  TGraph splineCurve;
  splineCurve.SetName("crystal_ball_mpv_spline_curve");
  splineCurve.SetLineColor(kGreen + 2);
  splineCurve.SetLineWidth(3);
  double splineMinimum = 0.0;
  double splineMaximum = 0.0;
  for (int point = 0; point < mpvGraph.GetN(); ++point) {
    double x = 0.0;
    double y = 0.0;
    mpvGraph.GetPoint(point, x, y);
    if (x >= kSplineMinimumTrackY && x <= kSplineMaximumTrackY) {
      splinePoints.SetPoint(splinePoints.GetN(), x, y);
    }
  }
  const bool hasSplineRange = splinePoints.GetN() >= 2;
  if (hasSplineRange) {
    mpvSpline = std::make_unique<TSpline3>(
        "crystal_ball_mpv_spline", &splinePoints, "b2e2");
    constexpr int kSplineSamples = 300;
    double endpointMpv = 0.0;
    splinePoints.GetPoint(0, splineMinimum, endpointMpv);
    splinePoints.GetPoint(splinePoints.GetN() - 1, splineMaximum, endpointMpv);
    for (int sample = 0; sample < kSplineSamples; ++sample) {
      const double x =
          splineMinimum +
          (splineMaximum - splineMinimum) * sample /
              (kSplineSamples - 1);
      splineCurve.SetPoint(sample, x, mpvSpline->Eval(x));
    }
    splineCurve.Draw("L SAME");
  } else {
    std::cerr << "Cannot build the MPV spline in trackY range ["
              << kSplineMinimumTrackY << ", " << kSplineMaximumTrackY
              << "]: at least two successful Crystal Ball fits in this range "
                 "are required."
              << std::endl;
  }
  mpvGraph.Draw("P SAME");
  TLegend legend(0.13, 0.82, 0.42, 0.89);
  legend.AddEntry(&mpvGraph, "Crystal Ball MPV", "p");
  if (mpvSpline) {
    legend.AddEntry(&splineCurve, "MPV spline fit", "l");
  }
  legend.Draw();
  if (mpvGraph.Write() <= 0 ||
      (mpvSpline &&
       (mpvSpline->Write() <= 0 || splineCurve.Write() <= 0)) ||
      canvas.Write() <= 0) {
    std::cerr << "Failed to write spatial-resolution output: " << rootPath
              << std::endl;
    return;
  }
  canvas.SaveAs(plotPath.c_str());

  if (!mpvSpline) {
    std::cerr << "Cannot calculate residuals without a valid MPV spline."
              << std::endl;
    return;
  }

  const std::string eventPath = runBase + ".root";
  std::unique_ptr<TFile> eventFile(TFile::Open(eventPath.c_str(), "READ"));
  if (!eventFile || eventFile->IsZombie()) {
    std::cerr << "Cannot open event data file: " << eventPath << std::endl;
    return;
  }
  auto* eventTree = dynamic_cast<TTree*>(eventFile->Get("results"));
  if (!eventTree) {
    std::cerr << "Tree 'results' not found in " << eventPath << std::endl;
    return;
  }

  double riseTime = 0.0;
  double trackY = 0.0;
  if (eventTree->SetBranchAddress("t_rise", &riseTime) < 0 ||
      eventTree->SetBranchAddress("y1", &trackY) < 0) {
    std::cerr << "Failed to bind 't_rise' and 'y1' in tree 'results'."
              << std::endl;
    eventTree->ResetBranchAddresses();
    return;
  }

  const double residualBinWidth =
      (kSplineMaximumTrackY - kSplineMinimumTrackY) / kResidualBins;
  std::array<std::unique_ptr<TH1D>, kResidualBins> residualHistograms;
  for (int bin = 0; bin < kResidualBins; ++bin) {
    const double binMinimum =
        kSplineMinimumTrackY + bin * residualBinWidth;
    const double binMaximum = binMinimum + residualBinWidth;
    std::ostringstream title;
    title << "Residuals for true track y in [" << std::fixed
          << std::setprecision(3) << binMinimum << ", " << binMaximum
          << (bin + 1 == kResidualBins ? "]" : ")");
    residualHistograms[bin] = std::make_unique<TH1D>(
        ("trackY_residual_bin_" + twoDigits(bin + 1)).c_str(),
        (title.str() + ";Reconstructed y - true y [cm];Events").c_str(),
        390, -0.65, 0.65);
    residualHistograms[bin]->SetDirectory(nullptr);
    residualHistograms[bin]->SetLineColor(kBlue + 1);
    residualHistograms[bin]->SetFillColor(kAzure - 9);
  }

  Long64_t reconstructedEntries = 0;
  Long64_t rejectedEntries = 0;
  for (Long64_t entry = 0; entry < eventTree->GetEntries(); ++entry) {
    if (eventTree->GetEntry(entry) <= 0) {
      std::cerr << "Failed to read event " << entry << " from " << eventPath
                << std::endl;
      eventTree->ResetBranchAddresses();
      return;
    }
    if (!std::isfinite(riseTime) || !std::isfinite(trackY) || riseTime < 0.0 ||
        trackY < kSplineMinimumTrackY ||
        trackY > kSplineMaximumTrackY) {
      continue;
    }

    double reconstructedTrackY = 0.0;
    if (!invertSpline(*mpvSpline, splineMinimum, splineMaximum, riseTime,
                      reconstructedTrackY)) {
      ++rejectedEntries;
      continue;
    }

    int bin = static_cast<int>(
        (trackY - kSplineMinimumTrackY) / residualBinWidth);
    if (bin == kResidualBins) {
      bin = kResidualBins - 1;
    }
    residualHistograms[bin]->Fill(reconstructedTrackY - trackY);
    ++reconstructedEntries;
  }
  eventTree->ResetBranchAddresses();

  TCanvas residualCanvas("c_trackY_residuals",
                         "Track y residuals by true position", 1500, 900);
  residualCanvas.Divide(4, 2);
  TGraphErrors residualSigmaGraph;
  residualSigmaGraph.SetName("residual_crystal_ball_sigma_vs_slice");
  residualSigmaGraph.SetTitle(
      "Residual width vs true track y;True track y slice center [cm];"
      "Width [#mum]");
  residualSigmaGraph.SetMarkerStyle(20);
  residualSigmaGraph.SetMarkerColor(kBlue + 1);
  residualSigmaGraph.SetLineColor(kBlue + 1);
  TGraph residualRmsGraph;
  residualRmsGraph.SetName("residual_rms_vs_slice");
  residualRmsGraph.SetTitle(
      "Residual RMS vs true track y;True track y slice center [cm];"
      "Width [#mum]");
  residualRmsGraph.SetMarkerStyle(21);
  residualRmsGraph.SetMarkerColor(kRed + 1);
  residualRmsGraph.SetLineColor(kRed + 1);
  std::array<std::unique_ptr<TLegend>, kResidualBins> residualLegends;
  for (int bin = 0; bin < kResidualBins; ++bin) {
    auto& residualHistogram = *residualHistograms[bin];
    const double sliceCenter =
        kSplineMinimumTrackY + (bin + 0.5) * residualBinWidth;
    const double residualRmsUm = residualHistogram.GetRMS() * 10000.0;
    residualRmsGraph.SetPoint(bin, sliceCenter, residualRmsUm);

    residualLegends[bin] =
        std::make_unique<TLegend>(0.52, 0.77, 0.89, 0.89);
    if (residualHistogram.GetEntries() > 0.0 &&
        residualHistogram.GetMaximum() > 0.0) {
      const double binWidth = residualHistogram.GetXaxis()->GetBinWidth(1);
      const double sigmaStart =
          std::max(residualHistogram.GetRMS(), binWidth);
      const double meanStart = residualHistogram.GetXaxis()->GetBinCenter(
          residualHistogram.GetMaximumBin());
      TF1 residualFit(("residual_crystal_ball_fit_" + twoDigits(bin + 1)).c_str(),
                      crystalBall, -0.3, 0.3, 5);
      residualFit.SetParNames("Amplitude", "MPV", "Sigma", "Alpha", "N");
      residualFit.SetParameters(residualHistogram.GetMaximum(), meanStart,
                                sigmaStart, 1.5, 3.0);
      residualFit.SetParLimits(0, 0.0,
                               residualHistogram.GetMaximum() * 100.0);
      residualFit.SetParLimits(1, -0.3, 0.3);
      residualFit.SetParLimits(2, binWidth * 0.1, 1.3);
      residualFit.SetParLimits(3, 0.1, 10.0);
      residualFit.SetParLimits(4, 1.01, 100.0);
      const int fitStatus = residualHistogram.Fit(&residualFit, "QRS0");
      if (fitStatus == 0) {
        auto* fittedFunction =
            residualHistogram.GetFunction(residualFit.GetName());
        if (fittedFunction) {
          fittedFunction->SetLineColor(kRed + 1);
          fittedFunction->SetLineWidth(3);
          residualLegends[bin]->AddEntry(fittedFunction, "Crystal Ball fit",
                                         "l");
          residualSigmaGraph.SetPoint(residualSigmaGraph.GetN(), sliceCenter,
                                      residualFit.GetParameter(2) * 10000.0);
          residualSigmaGraph.SetPointError(
              residualSigmaGraph.GetN() - 1, 0.0,
              residualFit.GetParError(2) * 10000.0);
        } else {
          std::cerr << "Crystal Ball function missing for residual slice "
                    << bin + 1 << "." << std::endl;
        }
      } else {
        std::cerr << "Crystal Ball fit failed for residual slice " << bin + 1
                  << " (status " << fitStatus << ")." << std::endl;
      }
    } else {
      std::cerr << "Residual slice " << bin + 1
                << " is empty; skipping Crystal Ball fit." << std::endl;
    }

    outputFile->cd();
    if (residualHistogram.Write() <= 0) {
      std::cerr << "Failed to write residual histogram for trackY bin "
                << bin + 1 << "." << std::endl;
      return;
    }
    residualCanvas.cd(bin + 1);
    residualHistogram.Draw("HIST");
    if (residualHistogram.GetFunction(
            ("residual_crystal_ball_fit_" + twoDigits(bin + 1)).c_str())) {
      auto* fittedFunction = residualHistogram.GetFunction(
          ("residual_crystal_ball_fit_" + twoDigits(bin + 1)).c_str());
      fittedFunction->Draw("SAME");
      residualLegends[bin]->Draw();
    }
  }
  outputFile->cd();
  if (residualCanvas.Write() <= 0) {
    std::cerr << "Failed to write residual-distribution canvas to " << rootPath
              << std::endl;
    return;
  }
  const std::string residualPlotPath =
      runBase + "_spatial_residuals.png";
  residualCanvas.SaveAs(residualPlotPath.c_str());

  residualRmsGraph.SetLineWidth(2);
  TMultiGraph residualWidths;
  residualWidths.SetName("residual_widths_vs_slice");
  residualWidths.SetTitle(
      "Residual widths by true track y;True track y slice center [cm];"
      "Width [#mum]");
  residualWidths.Add(&residualSigmaGraph, "P");
  residualWidths.Add(&residualRmsGraph, "LP");
  TCanvas widthsCanvas("c_residual_widths",
                       "Residual sigma and RMS by true position", 1000, 700);
  residualWidths.Draw("A");
  residualWidths.GetXaxis()->SetLimits(kSplineMinimumTrackY,
                                       kSplineMaximumTrackY);
  residualWidths.GetYaxis()->SetTitle("Width [#mum]");
  residualWidths.SetMinimum(0.0);
  residualWidths.SetMaximum(600.0);
  TLegend widthsLegend(0.62, 0.76, 0.88, 0.89);
  widthsLegend.AddEntry(&residualSigmaGraph, "Crystal Ball #sigma", "p");
  widthsLegend.AddEntry(&residualRmsGraph, "Residual RMS", "lp");
  widthsLegend.Draw();
  outputFile->cd();
  if (residualSigmaGraph.Write() <= 0 || residualRmsGraph.Write() <= 0 ||
      residualWidths.Write() <= 0 || widthsCanvas.Write() <= 0) {
    std::cerr << "Failed to write residual-width output to " << rootPath
              << std::endl;
    return;
  }
  const std::string widthsPlotPath =
      runBase + "_spatial_residual_widths.png";
  widthsCanvas.SaveAs(widthsPlotPath.c_str());
  outputFile->Close();

  std::ifstream plotFile(plotPath, std::ios::binary);
  if (!plotFile) {
    std::cerr << "Failed to create spatial-resolution plot: " << plotPath
              << std::endl;
    return;
  }
  std::ifstream timeResolutionPlotFile(timeResolutionPlotPath,
                                       std::ios::binary);
  if (!timeResolutionPlotFile) {
    std::cerr << "Failed to create time-resolution plot: "
              << timeResolutionPlotPath << std::endl;
    return;
  }
  std::ifstream residualPlotFile(residualPlotPath, std::ios::binary);
  if (!residualPlotFile) {
    std::cerr << "Failed to create residual-distribution plot: "
              << residualPlotPath << std::endl;
    return;
  }
  std::ifstream widthsPlotFile(widthsPlotPath, std::ios::binary);
  if (!widthsPlotFile) {
    std::cerr << "Failed to create residual-width plot: " << widthsPlotPath
              << std::endl;
    return;
  }
  std::cout << "Spatial-resolution plot saved to " << plotPath << std::endl;
  std::cout << "Crystal Ball time-resolution plot saved to "
            << timeResolutionPlotPath << " (" << timeResolutionGraph.GetN()
            << " successful slice fits)." << std::endl;
  std::cout << "MPV graph saved to " << rootPath << std::endl;
  std::cout << "Residual distributions saved to " << residualPlotPath
            << " (" << reconstructedEntries << " events; " << rejectedEntries
            << " events could not be inverted uniquely)." << std::endl;
  std::cout << "Residual sigma and RMS plot saved to " << widthsPlotPath
            << std::endl;
}
