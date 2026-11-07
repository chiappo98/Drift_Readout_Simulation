#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TTree.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

void DT_investigations(const char* runNumber, double tRiseMin,
                       double tRiseMax, double trackYMin, double trackYMax) {
  if (!runNumber || std::string(runNumber).empty()) {
    std::cerr << "A run number is required." << std::endl;
    return;
  }
  if (!std::isfinite(tRiseMin) || !std::isfinite(tRiseMax) ||
      !std::isfinite(trackYMin) || !std::isfinite(trackYMax) ||
      tRiseMin > tRiseMax || trackYMin > trackYMax) {
    std::cerr << "Selection limits must be finite and ordered." << std::endl;
    return;
  }

  const std::string runPath = std::string("results/run") + runNumber;
  const std::string resultsPath = runPath + ".root";
  const std::string waveformPath = runPath + "_waveform.root";
  const std::string outputPath = runPath + "_selectedWF.pdf";
  constexpr std::size_t maxWaveforms = 40;

  TFile resultsFile(resultsPath.c_str(), "READ");
  if (resultsFile.IsZombie()) {
    std::cerr << "Cannot open results file: " << resultsPath << std::endl;
    return;
  }
  auto* resultsTree = dynamic_cast<TTree*>(resultsFile.Get("results"));
  if (!resultsTree) {
    std::cerr << "Tree 'results' not found in " << resultsPath << std::endl;
    return;
  }

  TFile waveformFile(waveformPath.c_str(), "READ");
  if (waveformFile.IsZombie()) {
    std::cerr << "Cannot open waveform file: " << waveformPath << std::endl;
    return;
  }
  auto* waveformTree = dynamic_cast<TTree*>(waveformFile.Get("waveforms"));
  if (!waveformTree) {
    std::cerr << "Tree 'waveforms' not found in " << waveformPath << std::endl;
    return;
  }
  for (const char* branch : {"t_rise", "y1"}) {
    if (!resultsTree->GetBranch(branch)) {
      std::cerr << "Branch '" << branch << "' not found in " << resultsPath
                << std::endl;
      return;
    }
  }
  for (const char* branch : {"waveform", "raw_signal"}) {
    if (!waveformTree->GetBranch(branch)) {
      std::cerr << "Branch '" << branch << "' not found in " << waveformPath
                << ". Rerun the simulation to save both signals."
                << std::endl;
      return;
    }
  }

  const bool resultsHasEvent = resultsTree->GetBranch("event") != nullptr;
  const bool waveformsHaveEvent = waveformTree->GetBranch("event") != nullptr;
  const bool matchByEvent = resultsHasEvent && waveformsHaveEvent;
  if (!matchByEvent &&
      resultsTree->GetEntries() != waveformTree->GetEntries()) {
    std::cerr << "Cannot match entries by position: the results and waveform "
                 "trees have different entry counts."
              << std::endl;
    return;
  }

  double tRise = 0.0;
  double trackY = 0.0;
  ULong64_t resultsEvent = 0;
  if (resultsTree->SetBranchAddress("t_rise", &tRise) < 0 ||
      resultsTree->SetBranchAddress("y1", &trackY) < 0 ||
      (matchByEvent &&
       resultsTree->SetBranchAddress("event", &resultsEvent) < 0)) {
    std::cerr << "Failed to bind branches in " << resultsPath << std::endl;
    resultsTree->ResetBranchAddresses();
    return;
  }

  struct SelectedEvent {
    double tRise;
    double trackY;
  };
  std::unordered_map<ULong64_t, SelectedEvent> selectedEvents;
  for (Long64_t entry = 0; entry < resultsTree->GetEntries(); ++entry) {
    if (resultsTree->GetEntry(entry) <= 0) {
      std::cerr << "Failed to read results entry " << entry << std::endl;
      resultsTree->ResetBranchAddresses();
      return;
    }
    if (tRise < tRiseMin || tRise > tRiseMax || trackY < trackYMin ||
        trackY > trackYMax) {
      continue;
    }
    const ULong64_t eventId =
        matchByEvent ? resultsEvent : static_cast<ULong64_t>(entry);
    if (!selectedEvents.emplace(eventId, SelectedEvent{tRise, trackY}).second) {
      std::cerr << "Duplicate event number " << eventId << " in "
                << resultsPath << std::endl;
      resultsTree->ResetBranchAddresses();
      return;
    }
  }
  resultsTree->ResetBranchAddresses();
  if (selectedEvents.empty()) {
    std::cout << "No events satisfy the requested selection." << std::endl;
    return;
  }

  ULong64_t waveformEvent = 0;
  std::vector<float>* waveform = nullptr;
  std::vector<float>* rawSignal = nullptr;
  if ((matchByEvent &&
       waveformTree->SetBranchAddress("event", &waveformEvent) < 0) ||
      waveformTree->SetBranchAddress("waveform", &waveform) < 0 ||
      waveformTree->SetBranchAddress("raw_signal", &rawSignal) < 0) {
    std::cerr << "Failed to bind branches in " << waveformPath << std::endl;
    waveformTree->ResetBranchAddresses();
    return;
  }

  const std::string canvasName =
      std::string("c_DT_investigations_run") + runNumber;
  auto* canvas = new TCanvas(canvasName.c_str(), "Selected waveforms", 1100, 900);
  std::size_t plottedEvents = 0;
  std::size_t matchedEvents = 0;
  for (Long64_t entry = 0; entry < waveformTree->GetEntries(); ++entry) {
    if (waveformTree->GetEntry(entry) <= 0) {
      std::cerr << "Failed to read waveform entry " << entry << std::endl;
      waveformTree->ResetBranchAddresses();
      return;
    }
    const ULong64_t eventId =
        matchByEvent ? waveformEvent : static_cast<ULong64_t>(entry);
    const auto selected = selectedEvents.find(eventId);
    if (selected == selectedEvents.end()) {
      continue;
    }
    if (!waveform || waveform->empty()) {
      std::cerr << "Waveform for event " << eventId << " is empty."
                << std::endl;
      waveformTree->ResetBranchAddresses();
      return;
    }
    if (!rawSignal || rawSignal->empty() ||
        rawSignal->size() != waveform->size()) {
      std::cerr << "Raw signal for event " << eventId
                << " is empty or has a different number of samples."
                << std::endl;
      waveformTree->ResetBranchAddresses();
      return;
    }
    ++matchedEvents;
    if (plottedEvents >= maxWaveforms) {
      continue;
    }

    std::vector<double> times(waveform->size());
    std::vector<double> shapedSamples(waveform->size());
    std::vector<double> rawSamples(rawSignal->size());
    for (std::size_t sample = 0; sample < waveform->size(); ++sample) {
      times[sample] = static_cast<double>(sample);
      shapedSamples[sample] = (*waveform)[sample];
      rawSamples[sample] = (*rawSignal)[sample];
    }
    const auto& selectedEvent = selected->second;
    std::ostringstream shapedTitle;
    shapedTitle << "Event " << eventId << ", trackY = " << std::fixed
                << std::setprecision(2) << selectedEvent.trackY
                << " cm;Time [ns];Shaped signal [native units]";
    TGraph shapedGraph(static_cast<int>(waveform->size()), times.data(),
                       shapedSamples.data());
    shapedGraph.SetTitle(shapedTitle.str().c_str());
    shapedGraph.GetXaxis()->SetLimits(0.0, 1000.0);
    shapedGraph.SetLineColor(kBlue + static_cast<int>(plottedEvents % 5));
    shapedGraph.SetLineWidth(1);
    canvas->Clear();
    canvas->Divide(1, 2);
    TPad* shapedPad = static_cast<TPad*>(canvas->cd(1));
    shapedGraph.Draw("AL");
    shapedPad->Update();
    TLine tRiseLine(selectedEvent.tRise, shapedPad->GetUymin(),
                    selectedEvent.tRise, shapedPad->GetUymax());
    tRiseLine.SetLineColor(kRed);
    tRiseLine.SetLineStyle(2);
    tRiseLine.SetLineWidth(2);
    tRiseLine.Draw("same");

    auto* riseText = new TPaveText(0.55, 0.62, 0.85, 0.82, "NDC");
    riseText->SetBorderSize(0);
    riseText->SetFillColor(0);
    riseText->SetTextAlign(12);
    riseText->AddText(Form("t_{rise} = %.2f ns", selectedEvent.tRise));
    riseText->AddText(Form("distance = %.2f cm", selectedEvent.trackY));
    riseText->Draw();

    TGraph rawGraph(static_cast<int>(rawSignal->size()), times.data(),
                    rawSamples.data());
    rawGraph.SetTitle("Raw signal;Time [ns];Current [native units]");
    rawGraph.GetXaxis()->SetLimits(0.0, 1000.0);
    rawGraph.SetLineColor(kRed + 1);
    rawGraph.SetLineWidth(1);
    canvas->cd(2);
    rawGraph.Draw("AL");

    if (plottedEvents == 0) {
      canvas->Print((outputPath + "[").c_str());
    }
    canvas->Print(outputPath.c_str());
    ++plottedEvents;
  }
  waveformTree->ResetBranchAddresses();

  if (plottedEvents > 0) {
    canvas->Print((outputPath + "]").c_str());
  }
  if (plottedEvents == 0) {
    std::cerr << "No selected events were found in " << waveformPath
              << std::endl;
    return;
  }

  std::cout << "Plotted " << plottedEvents << " of " << matchedEvents
            << " matching events (maximum " << maxWaveforms << ")."
            << std::endl;
  if (matchedEvents != selectedEvents.size()) {
    std::cerr << "Some selected events have no matching waveform entry."
              << std::endl;
  }
  std::cout << "Waveform plot saved to " << outputPath << std::endl;
}
