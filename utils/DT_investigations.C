#include <TCanvas.h>
#include <TBox.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TEllipse.h>
#include <TH1D.h>
#include <TLine.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveText.h>
#include <TH2D.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
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
  const double absTrackYMin =
      trackYMin <= 0.0 && trackYMax >= 0.0
          ? 0.0
          : std::min(std::abs(trackYMin), std::abs(trackYMax));
  const double absTrackYMax =
      std::max(std::abs(trackYMin), std::abs(trackYMax));
  std::ostringstream outputName;
  outputName << runPath << "_selectedWF_" << std::setprecision(6)
             << absTrackYMin << "-" << absTrackYMax << "_" << tRiseMin << "-"
             << tRiseMax << ".pdf";
  const std::string outputPath = outputName.str();
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
  for (const char* branch : {"totalCharge", "t_rise", "x1", "y1", "z1",
                             "dx1", "dy1", "dz1", "clusterX", "clusterY",
                             "clusterZ", "clusterElectrons"}) {
    if (!resultsTree->GetBranch(branch)) {
      std::cerr << "Branch '" << branch << "' not found in " << resultsPath
                << std::endl;
      return;
    }
  }
  const auto bindOptionalVector = [&resultsTree](const char* primary,
                                                 const char* fallback,
                                                 std::vector<double>** branchData) {
    if (resultsTree->GetBranch(primary)) {
      return resultsTree->SetBranchAddress(primary, branchData) >= 0;
    }
    if (resultsTree->GetBranch(fallback)) {
      return resultsTree->SetBranchAddress(fallback, branchData) >= 0;
    }
    return false;
  };
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

  double totalCharge = 0.0;
  double tRise = 0.0;
  double trackX = 0.0;
  double trackY = 0.0;
  double trackZ = 0.0;
  double directionX = 0.0;
  double directionY = 0.0;
  double directionZ = 0.0;
  std::vector<double>* clusterX = nullptr;
  std::vector<double>* clusterY = nullptr;
  std::vector<double>* clusterZ = nullptr;
  std::vector<int>* clusterElectrons = nullptr;
  std::vector<double>* electronX = nullptr;
  std::vector<double>* electronY = nullptr;
  const bool hasElectronXY = bindOptionalVector("electron_x", "clusterElectronX",
                                               &electronX) &&
      bindOptionalVector("electron_y", "clusterElectronY", &electronY);
  ULong64_t resultsEvent = 0;
  if (resultsTree->SetBranchAddress("totalCharge", &totalCharge) < 0 ||
      resultsTree->SetBranchAddress("t_rise", &tRise) < 0 ||
      resultsTree->SetBranchAddress("x1", &trackX) < 0 ||
      resultsTree->SetBranchAddress("y1", &trackY) < 0 ||
      resultsTree->SetBranchAddress("z1", &trackZ) < 0 ||
      resultsTree->SetBranchAddress("dx1", &directionX) < 0 ||
      resultsTree->SetBranchAddress("dy1", &directionY) < 0 ||
      resultsTree->SetBranchAddress("dz1", &directionZ) < 0 ||
      resultsTree->SetBranchAddress("clusterX", &clusterX) < 0 ||
      resultsTree->SetBranchAddress("clusterY", &clusterY) < 0 ||
      resultsTree->SetBranchAddress("clusterZ", &clusterZ) < 0 ||
      resultsTree->SetBranchAddress("clusterElectrons", &clusterElectrons) < 0 ||
      (matchByEvent &&
       resultsTree->SetBranchAddress("event", &resultsEvent) < 0) ||
      (!hasElectronXY &&
       (resultsTree->GetBranch("electron_x") ||
        resultsTree->GetBranch("clusterElectronX") ||
        resultsTree->GetBranch("electron_y") ||
        resultsTree->GetBranch("clusterElectronY")))) {
    std::cerr << "Failed to bind branches in " << resultsPath << std::endl;
    resultsTree->ResetBranchAddresses();
    return;
  }

  struct SelectedEvent {
    double totalCharge;
    double tRise;
    double trackY;
    double trackX;
    double trackZ;
    double directionX;
    double directionY;
    double directionZ;
    std::vector<double> clusterX;
    std::vector<double> clusterY;
    std::vector<double> clusterZ;
    std::vector<int> clusterElectrons;
    std::vector<double> electronX;
    std::vector<double> electronY;
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
    if (!clusterX || !clusterY || !clusterZ || !clusterElectrons ||
        clusterX->size() != clusterY->size() ||
        clusterX->size() != clusterZ->size() ||
        clusterX->size() != clusterElectrons->size()) {
      std::cerr << "Cluster branches are empty or have mismatched sizes at "
                   "results entry "
                << entry << std::endl;
      resultsTree->ResetBranchAddresses();
      return;
    }
    const ULong64_t eventId =
        matchByEvent ? resultsEvent : static_cast<ULong64_t>(entry);
    std::vector<double> eventElectronX;
    std::vector<double> eventElectronY;
    if (electronX && electronY) {
      eventElectronX = *electronX;
      eventElectronY = *electronY;
    }
    SelectedEvent event{totalCharge, tRise, trackY, trackX, trackZ, directionX,
                        directionY, directionZ, *clusterX, *clusterY,
                        *clusterZ, *clusterElectrons, eventElectronX,
                        eventElectronY};
    if (!selectedEvents.emplace(eventId, std::move(event)).second) {
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
  canvas->Print((outputPath + "[").c_str());
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
    const double expectedTRise = std::abs(selectedEvent.trackY) / 0.005;
    TBox expectedTRiseBand(expectedTRise - 5.0, shapedPad->GetUymin(),
                           expectedTRise + 5.0, shapedPad->GetUymax());
    expectedTRiseBand.SetFillColorAlpha(kOrange + 1, 0.3);
    expectedTRiseBand.SetLineColor(kOrange + 1);
    expectedTRiseBand.SetLineWidth(0);
    expectedTRiseBand.Draw("same");
    TLine tRiseLine(selectedEvent.tRise, shapedPad->GetUymin(),
                    selectedEvent.tRise, shapedPad->GetUymax());
    tRiseLine.SetLineColor(kRed);
    tRiseLine.SetLineStyle(2);
    tRiseLine.SetLineWidth(2);
    tRiseLine.Draw("same");

    auto* riseLegend = new TLegend(0.57, 0.7, 0.87, 0.82);
    riseLegend->SetBorderSize(0);
    riseLegend->SetFillStyle(0);
    riseLegend->AddEntry(&expectedTRiseBand, "expected t_{rise}", "f");
    riseLegend->AddEntry(&tRiseLine, "measured t_{rise}", "l");
    riseLegend->Draw();

    auto* riseText = new TPaveText(0.55, 0.52, 0.85, 0.68, "NDC");
    riseText->SetBorderSize(0);
    riseText->SetFillColor(0);
    riseText->SetTextAlign(12);
    riseText->AddText(
        Form("measured t_{rise} = %.2f ns", selectedEvent.tRise));
    riseText->AddText(
        Form("distance = %.2f cm", std::abs(selectedEvent.trackY)));
    riseText->Draw();

    TGraph rawGraph(static_cast<int>(rawSignal->size()), times.data(),
                    rawSamples.data());
    rawGraph.SetTitle("Raw signal;Time [ns];Current [native units]");
    rawGraph.GetXaxis()->SetLimits(0.0, 1000.0);
    rawGraph.SetLineColor(kRed + 1);
    rawGraph.SetLineWidth(1);
    canvas->cd(2);
    rawGraph.Draw("AL");

    canvas->Print(outputPath.c_str());

    canvas->Clear();
    canvas->cd();
    const double directionMagnitude =
        std::sqrt(selectedEvent.directionX * selectedEvent.directionX +
                  selectedEvent.directionY * selectedEvent.directionY +
                  selectedEvent.directionZ * selectedEvent.directionZ);
    double trackEndX = selectedEvent.trackX;
    double trackEndY = selectedEvent.trackY;
    double maxTrackDistance = 0.0;
    if (directionMagnitude > 0.0) {
      for (std::size_t cluster = 0;
           cluster < selectedEvent.clusterX.size(); ++cluster) {
        const double dx = selectedEvent.clusterX[cluster] - selectedEvent.trackX;
        const double dy = selectedEvent.clusterY[cluster] - selectedEvent.trackY;
        const double dz = selectedEvent.clusterZ[cluster] - selectedEvent.trackZ;
        maxTrackDistance =
            std::max(maxTrackDistance,
                     (dx * selectedEvent.directionX +
                      dy * selectedEvent.directionY +
                      dz * selectedEvent.directionZ) /
                         directionMagnitude);
      }
      trackEndX += maxTrackDistance * selectedEvent.directionX /
                   directionMagnitude;
      trackEndY += maxTrackDistance * selectedEvent.directionY /
                   directionMagnitude;
    }

    std::ostringstream trackTitle;
    trackTitle << "Event " << eventId << " track and clusters;x [cm];y [cm]";
    const std::string frameName = "trackFrame_" + std::to_string(eventId);

    const bool hasElectronData = !selectedEvent.electronX.empty() &&
                                 !selectedEvent.electronY.empty() &&
                                 selectedEvent.electronX.size() ==
                                     selectedEvent.electronY.size();
    std::vector<double> fillX;
    std::vector<double> fillY;
    if (hasElectronData) {
      fillX = selectedEvent.electronX;
      fillY = selectedEvent.electronY;
    }
    if ((!hasElectronData || fillX.empty() || fillY.empty()) &&
        !selectedEvent.clusterX.empty() && !selectedEvent.clusterY.empty()) {
      fillX = selectedEvent.clusterX;
      fillY = selectedEvent.clusterY;
    }

    double frameXMin = -0.6;
    double frameXMax = 0.6;
    double frameYMin = -0.7;
    double frameYMax = 0.7;
    std::size_t validPoints = 0;
    if (fillX.size() == fillY.size()) {
      for (std::size_t i = 0; i < fillX.size() && i < fillY.size(); ++i) {
        if (!std::isfinite(fillX[i]) || !std::isfinite(fillY[i])) {
          continue;
        }
        frameXMin = std::min(frameXMin, fillX[i]);
        frameXMax = std::max(frameXMax, fillX[i]);
        frameYMin = std::min(frameYMin, fillY[i]);
        frameYMax = std::max(frameYMax, fillY[i]);
        ++validPoints;
      }
      if (validPoints > 0) {
        const double xPadMargin = std::max(0.05, 0.10 * (frameXMax - frameXMin));
        const double yPadMargin = std::max(0.05, 0.10 * (frameYMax - frameYMin));
        frameXMin -= xPadMargin;
        frameXMax += xPadMargin;
        frameYMin -= yPadMargin;
        frameYMax += yPadMargin;
      }
    }

    TH2D trackFrame(frameName.c_str(), trackTitle.str().c_str(), 100,
                    frameXMin, frameXMax, 100, frameYMin, frameYMax);
    trackFrame.SetDirectory(nullptr);
    trackFrame.SetStats(false);

    canvas->Clear();
    canvas->cd();
    canvas->SetLeftMargin(0.10);
    canvas->SetRightMargin(0.10);
    canvas->SetTopMargin(0.08);
    canvas->SetBottomMargin(0.12);
    trackFrame.Draw("AXIS");

    std::unique_ptr<TH2D> electronBackground;
    if (fillX.size() != fillY.size()) {
      std::cerr << "Event " << eventId
                << " has mismatched electron/cluster coordinate arrays; "
                   "skipping the electron distribution plot."
                << std::endl;
      auto* noDataText = new TPaveText(0.15, 0.35, 0.85, 0.65, "NDC");
      noDataText->SetBorderSize(0);
      noDataText->SetFillStyle(0);
      noDataText->SetTextAlign(22);
      noDataText->AddText("No valid electron positions available");
      noDataText->Draw();
    } else if (validPoints == 0) {
      std::cerr << "Event " << eventId
                << " has zero valid finite coordinates for the electron "
                   "distribution plot."
                << std::endl;
      auto* noDataText = new TPaveText(0.15, 0.35, 0.85, 0.65, "NDC");
      noDataText->SetBorderSize(0);
      noDataText->SetFillStyle(0);
      noDataText->SetTextAlign(22);
      noDataText->AddText("No valid electron positions after filtering");
      noDataText->Draw();
    } else {
      std::cout << "Event " << eventId
                << ": electron points valid=" << validPoints
                << ", x in [" << frameXMin << ", " << frameXMax << "] cm, y in ["
                << frameYMin << ", " << frameYMax << "] cm" << std::endl;
      electronBackground = std::make_unique<TH2D>(
          (frameName + "_electrons").c_str(),
          Form("Event %llu electron distribution; x [cm]; y [cm]", eventId),
          100, frameXMin, frameXMax, 100, frameYMin, frameYMax);
      electronBackground->SetDirectory(nullptr);
      electronBackground->SetStats(false);
      for (std::size_t i = 0; i < fillX.size() && i < fillY.size(); ++i) {
        if (!std::isfinite(fillX[i]) || !std::isfinite(fillY[i])) {
          continue;
        }
        electronBackground->Fill(fillX[i], fillY[i]);
      }
      const double maxBinContent = std::max(1.0, electronBackground->GetMaximum());
      electronBackground->SetMinimum(0.0);
      electronBackground->SetMaximum(maxBinContent);
      electronBackground->GetZaxis()->SetRangeUser(0.0, maxBinContent);
      electronBackground->SetContour(50);
      electronBackground->Draw("COLZ");
      canvas->Modified();
      canvas->Update();
      auto* palette = static_cast<TPaletteAxis*>(
          electronBackground->GetListOfFunctions()->FindObject("palette"));
      if (palette) {
        palette->SetX1NDC(0.88);
        palette->SetX2NDC(0.93);
        palette->SetY1NDC(0.12);
        palette->SetY2NDC(0.90);
        palette->SetTitle("electrons/bin");
      }
    }

    canvas->cd();
    TLine leftBoundary(-0.6, frameYMin, -0.6, frameYMax);
    leftBoundary.SetLineColor(kBlack);
    leftBoundary.SetLineStyle(2);
    leftBoundary.SetLineWidth(1);
    leftBoundary.Draw("same");

    TLine rightBoundary(0.6, frameYMin, 0.6, frameYMax);
    rightBoundary.SetLineColor(kBlack);
    rightBoundary.SetLineStyle(2);
    rightBoundary.SetLineWidth(1);
    rightBoundary.Draw("same");

    TLine trackLine(selectedEvent.trackX, selectedEvent.trackY, trackEndX,
                    trackEndY);
    trackLine.SetLineColor(kBlue + 1);
    trackLine.SetLineWidth(2);
    if (maxTrackDistance > 0.0) {
      trackLine.Draw("same");
    }
    TMarker trackStart(selectedEvent.trackX, selectedEvent.trackY, 29);
    trackStart.SetMarkerColor(kBlack);
    trackStart.SetMarkerSize(1.7);
    trackStart.Draw("same");

    const int maxElectrons =
        selectedEvent.clusterElectrons.empty()
            ? 0
            : *std::max_element(selectedEvent.clusterElectrons.begin(),
                                selectedEvent.clusterElectrons.end());
    long long totalPrimaryElectrons = 0;
    for (const int electrons : selectedEvent.clusterElectrons) {
      totalPrimaryElectrons += electrons;
    }
    std::vector<std::unique_ptr<TMarker>> clusterMarkers;
    clusterMarkers.reserve(selectedEvent.clusterX.size());
    for (std::size_t cluster = 0;
         cluster < selectedEvent.clusterX.size(); ++cluster) {
      const int electrons = selectedEvent.clusterElectrons[cluster];
      if (electrons <= 0 || maxElectrons <= 0) {
        continue;
      }
      auto marker = std::make_unique<TMarker>(
          selectedEvent.clusterX[cluster], selectedEvent.clusterY[cluster],
          24);
      marker->SetMarkerColor(kRed + 1);
      marker->SetMarkerSize(
          5 * std::sqrt(static_cast<double>(electrons) / maxElectrons));
      marker->Draw("same");
      clusterMarkers.push_back(std::move(marker));
    }

    const std::array<double, 2> referenceY{{-0.7, 0.7}};
    std::vector<std::unique_ptr<TEllipse>> referenceCircles;
    referenceCircles.reserve(referenceY.size() + 1);
    for (const double y : referenceY) {
      auto circle = std::make_unique<TEllipse>(0.0, y, 0.0125, 0.0125);
      circle->SetLineColor(kBlack);
      circle->SetLineWidth(2);
      circle->SetFillStyle(0);
      circle->Draw("same");
      referenceCircles.push_back(std::move(circle));
    }
    auto originCircle = std::make_unique<TEllipse>(0.0, 0.0, 0.0075, 0.0075);
    originCircle->SetLineColor(kBlack);
    originCircle->SetFillColor(kBlack);
    originCircle->SetFillStyle(1001);
    originCircle->Draw("same");
    referenceCircles.push_back(std::move(originCircle));

    auto* trackInfo = new TPaveText(0.50, 0.6, 0.85, 0.8, "NDC");
    trackInfo->SetBorderSize(0);
    trackInfo->SetFillColor(0);
    trackInfo->SetTextAlign(12);
    trackInfo->AddText(Form("Start: (%.3f, %.3f, %.3f) cm",
                            selectedEvent.trackX, selectedEvent.trackY,
                            selectedEvent.trackZ));
    trackInfo->AddText(
        Form("Clusters generated: %llu",
             static_cast<unsigned long long>(
                 selectedEvent.clusterElectrons.size())));
    trackInfo->AddText(
        Form("Max electrons in a cluster: %d (%.3f fC)", maxElectrons,
             maxElectrons * 1.602e-4));
    trackInfo->AddText(
        Form("Total primary electrons: %lld (%.3f fC)", totalPrimaryElectrons,
             totalPrimaryElectrons * 1.602e-4));
    trackInfo->AddText(
        Form("Total charge: %.2f fC", selectedEvent.totalCharge));
    trackInfo->AddText("Cluster circle area #propto N_{e}");
    trackInfo->Draw();
    canvas->Modified();
    canvas->Update();
    canvas->Print(outputPath.c_str());
    ++plottedEvents;
  }
  waveformTree->ResetBranchAddresses();

  canvas->Clear();
  canvas->cd();
  TH1D electronRatioHistogram(
      ("h_electron_ratio_run" + std::string(runNumber)).c_str(),
      "Maximum cluster electrons / total primary electrons;Ratio (0 if total = 0);Selected events",
      100, 0.0, 1.001);
  electronRatioHistogram.SetDirectory(nullptr);
  electronRatioHistogram.SetFillColor(kBlue - 7);
  for (const auto& eventEntry : selectedEvents) {
    const auto& electrons = eventEntry.second.clusterElectrons;
    const int eventMaxElectrons =
        electrons.empty() ? 0
                          : *std::max_element(electrons.begin(), electrons.end());
    long long eventTotalElectrons = 0;
    for (const int count : electrons) {
      eventTotalElectrons += count;
    }
    const double ratio =
        eventTotalElectrons > 0
            ? static_cast<double>(eventMaxElectrons) / eventTotalElectrons
            : 0.0;
    electronRatioHistogram.Fill(ratio);
  }
  electronRatioHistogram.Draw();
  canvas->Modified();
  canvas->Update();
  canvas->Print(outputPath.c_str());
  canvas->Print((outputPath + "]").c_str());

  if (plottedEvents == 0) {
    std::cerr << "No selected events were found in " << waveformPath
              << "; the PDF contains the ratio histogram for all selected "
                 "results events."
              << std::endl;
  }

  std::cout << "Plotted " << plottedEvents << " of " << matchedEvents
            << " matching events (maximum " << maxWaveforms << ")."
            << std::endl;
  if (matchedEvents != selectedEvents.size()) {
    std::cerr << "Some selected events have no matching waveform entry."
              << std::endl;
  }
  std::cout << "Investigation plots saved to " << outputPath << std::endl;
  
}
