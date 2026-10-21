#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TPad.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

void analyse_correlations(const char* runNumber) {
  if (!runNumber || std::string(runNumber).empty()) {
    std::cerr << "A run number is required for correlation analysis." << std::endl;
    return;
  }

  const std::string runPath = std::string("results/run") + runNumber;
  const std::string inputPath = runPath + ".root";
  const std::string outputPath = runPath + "_correlations.root";
  const std::string imagePath = runPath + "_correlations.png";

  TFile inputFile(inputPath.c_str(), "READ");
  if (inputFile.IsZombie()) {
    std::cerr << "Cannot open input file for correlations: " << inputPath << std::endl;
    return;
  }

  auto* tree = dynamic_cast<TTree*>(inputFile.Get("results"));
  if (!tree) {
    std::cerr << "Tree 'results' not found in " << inputPath << std::endl;
    return;
  }

  const std::array<const char*, 6> requiredBranches{{
      "totalCharge", "q_peak", "t_rise", "t_fall", "t_peak", "y1"}};
  for (const auto* branch : requiredBranches) {
    if (!tree->GetBranch(branch)) {
      std::cerr << "Branch '" << branch
                << "' not found for correlation analysis." << std::endl;
      return;
    }
  }

  double totalCharge = 0.0;
  double peakAmplitude = 0.0;
  double riseTime = 0.0;
  double fallTime = 0.0;
  double peakTime = 0.0;
  double trackY = 0.0;
  if (tree->SetBranchAddress("totalCharge", &totalCharge) < 0 ||
      tree->SetBranchAddress("q_peak", &peakAmplitude) < 0 ||
      tree->SetBranchAddress("t_rise", &riseTime) < 0 ||
      tree->SetBranchAddress("t_fall", &fallTime) < 0 ||
      tree->SetBranchAddress("t_peak", &peakTime) < 0 ||
      tree->SetBranchAddress("y1", &trackY) < 0) {
    std::cerr << "Failed to bind branches for correlation analysis." << std::endl;
    tree->ResetBranchAddresses();
    return;
  }

  const std::array<const char*, 10> names{{
      "h2_totalCharge_vs_peakAmplitude",
      "h2_totalCharge_vs_fallRise",
      "h2_peakAmplitude_vs_fallRise",
      "h2_totalCharge_vs_peakRise",
      "h2_peakAmplitude_vs_peakRise",
      "h2_peakAmplitudeOverTotalCharge_vs_trackY",
      "h2_fallRise_vs_trackY",
      "h2_peakRise_vs_trackY",
      "h2_riseTime_vs_peakTime",
      "h2_riseTime_vs_trackY"}};
  const std::array<const char*, 10> titles{{
      "Peak amplitude vs total charge;Total charge [fC];Peak amplitude [fC]",
      "ToT vs total charge;Total charge [fC];ToT [ns]",
      "ToT vs peak amplitude;Peak amplitude [fC];ToT [ns]",
      "Total charge vs peak time - time-at-threshold;Total charge [fC];t_{peak} - t_{thr} [ns]",
      "Peak amplitude vs peak time - time-at-threshold;Peak amplitude [fC];t_{peak} - t_{thr} [ns]",
      "Peak amplitude / total charge vs track y position;y [cm];Peak amplitude / total charge",
      "ToT vs track y position;y [cm];ToT [ns]",
      "Peak time - time-at-threshold vs track y position;y [cm];t_{peak} - t_{thr} [ns]",
      "Time-at-threshold vs peak time;Time-at-threshold [ns];t_{peak} [ns]",
      "Time-at-threshold vs track y position;y [cm];Time-at-threshold [ns]"}};
  std::array<std::vector<std::pair<double, double>>, 10> points;

  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    if (tree->GetEntry(entry) <= 0) {
      std::cerr << "Failed to read entry " << entry
                << " for correlation analysis." << std::endl;
      tree->ResetBranchAddresses();
      return;
    }

    const double amplitude = std::abs(peakAmplitude);
    points[0].emplace_back(totalCharge, amplitude);

    if (totalCharge != 0.0) {
      points[5].emplace_back(trackY, amplitude / totalCharge);
    }
    if (riseTime < 0.0) {
      continue;
    }

    points[8].emplace_back(riseTime, peakTime);
    points[9].emplace_back(trackY, riseTime);
    const double peakRise = peakTime - riseTime;
    points[3].emplace_back(totalCharge, peakRise);
    points[4].emplace_back(amplitude, peakRise);
    points[7].emplace_back(trackY, peakRise);

    if (fallTime < 0.0) {
      continue;
    }
    const double fallRise = fallTime - riseTime;
    points[1].emplace_back(totalCharge, fallRise);
    points[2].emplace_back(amplitude, fallRise);
    points[6].emplace_back(trackY, fallRise);
  }
  tree->ResetBranchAddresses();

  const auto axisRange = [](const std::vector<double>& values) {
    if (values.empty()) {
      return std::make_pair(0.0, 1.0);
    }
    const auto bounds = std::minmax_element(values.begin(), values.end());
    double minimum = *bounds.first;
    double maximum = *bounds.second;
    if (minimum == maximum) {
      const double padding = std::max(0.5, std::abs(minimum) * 0.05);
      return std::make_pair(minimum - padding, maximum + padding);
    }
    const double padding = (maximum - minimum) * 0.05;
    return std::make_pair(minimum - padding, maximum + padding);
  };

  TFile outputFile(outputPath.c_str(), "RECREATE");
  if (outputFile.IsZombie()) {
    std::cerr << "Cannot create correlation file: " << outputPath << std::endl;
    return;
  }
  outputFile.cd();
  TCanvas canvas("c_correlations", "Simulation correlations", 1800, 1000);
  canvas.Divide(5, 2);

  for (std::size_t i = 0; i < points.size(); ++i) {
    std::vector<double> xValues;
    std::vector<double> yValues;
    xValues.reserve(points[i].size());
    yValues.reserve(points[i].size());
    for (const auto& point : points[i]) {
      xValues.push_back(point.first);
      yValues.push_back(point.second);
    }
    const auto xRange = axisRange(xValues);
    const auto yRange = axisRange(yValues);
    auto* histogram = new TH2D(names[i], titles[i], 100, xRange.first,
                               xRange.second, 100, yRange.first,
                               yRange.second);
    for (const auto& point : points[i]) {
      histogram->Fill(point.first, point.second);
    }
    canvas.cd(static_cast<int>(i + 1));
    histogram->Draw("COLZ");
    histogram->Write();
  }

  canvas.Write();
  canvas.SaveAs(imagePath.c_str());
  std::cout << "Correlation histograms saved to " << outputPath << std::endl;
  std::cout << "Correlation canvas saved to " << imagePath << std::endl;
}

void analyse_results(const char* runNumber) {
  if (!runNumber || std::string(runNumber).empty()) {
    std::cerr << "A run number is required." << std::endl;
    return;
  }

  const std::string runPath = std::string("results/run") + runNumber;
  const std::string inputPath = runPath + ".root";
  const std::string outputPath = runPath + "_histograms.root";
  const std::string imagePath = runPath + "_histograms.png";

  TFile inputFile(inputPath.c_str(), "READ");
  if (inputFile.IsZombie()) {
    std::cerr << "Cannot open input file: " << inputPath << std::endl;
    return;
  }

  auto* tree = dynamic_cast<TTree*>(inputFile.Get("results"));
  if (!tree) {
    std::cerr << "Tree 'results' not found in " << inputPath << std::endl;
    return;
  }

  struct HistogramSpec {
    const char* expression;
    const char* name;
    const char* title;
  };
  const std::array<const char*, 17> requiredBranches{{
      "totalCharge", "q_peak", "t_peak", "t_rise", "t_fall",
      "x1", "y1", "z1", "dx1", "dy1", "dz1",
      "clusterX", "clusterY", "clusterZ", "clusterT", "clusterE",
      "clusterElectrons"}};
  const std::array<HistogramSpec, 14> specs{{
      {"totalCharge", "h_totalCharge", "Total charge;Charge [fC];Events"},
      {"abs(q_peak)", "h_q_peak", "Peak amplitude;Amplitude [fC];Events"},
      {"t_peak", "h_t_peak", "Peak time;Time [ns];Events"},
      {"t_rise", "h_t_rise", "Rise threshold crossing;Time [ns];Events"},
      {"t_fall", "h_t_fall", "Fall threshold crossing;Time [ns];Events"},
      {"x1", "h_x1", "Track start x;Position [cm];Events"},
      {"y1", "h_y1", "Track start y;Position [cm];Events"},
      {"z1", "h_z1", "Track start z;Position [cm];Events"},
      {"dx1", "h_dx1", "Track direction dx;Direction cosine;Events"},
      {"dy1", "h_dy1", "Track direction dy;Direction cosine;Events"},
      {"dz1", "h_dz1", "Track direction dz;Direction cosine;Events"},
      {"totalCharge != 0 ? abs(q_peak) / totalCharge : 0", "h_q_peak_over_totalCharge",
       "Peak/total charge ratio;Q peak / Q total;Events"},
      {"t_peak - t_rise", "h_t_peak_minus_t_rise",
       "Peak time - rise time;Time difference [ns];Events"},
      {"t_fall - t_rise", "h_t_fall_minus_t_rise",
       "Time-over-Threshold;ToT[ns];Events"}}};

  for (const auto* branch : requiredBranches) {
    if (!tree->GetBranch(branch)) {
      std::cerr << "Branch '" << branch << "' not found in tree 'results'." << std::endl;
      return;
    }
  }

  TFile outputFile(outputPath.c_str(), "RECREATE");
  if (outputFile.IsZombie()) {
    std::cerr << "Cannot create output file: " << outputPath << std::endl;
    return;
  }

  outputFile.cd();
  TCanvas canvas("c_results", "Simulation results", 1600, 1200);
  canvas.Divide(4, 5);

  double totalCharge = 0.0;
  double peakAmplitude = 0.0;
  double peakTime = 0.0;
  double riseTime = 0.0;
  double fallTime = 0.0;
  double trackX = 0.0;
  double trackY = 0.0;
  double trackZ = 0.0;
  double directionX = 0.0;
  double directionY = 0.0;
  double directionZ = 0.0;
  if (tree->SetBranchAddress("totalCharge", &totalCharge) < 0 ||
      tree->SetBranchAddress("q_peak", &peakAmplitude) < 0 ||
      tree->SetBranchAddress("t_peak", &peakTime) < 0 ||
      tree->SetBranchAddress("t_rise", &riseTime) < 0 ||
      tree->SetBranchAddress("t_fall", &fallTime) < 0 ||
      tree->SetBranchAddress("x1", &trackX) < 0 ||
      tree->SetBranchAddress("y1", &trackY) < 0 ||
      tree->SetBranchAddress("z1", &trackZ) < 0 ||
      tree->SetBranchAddress("dx1", &directionX) < 0 ||
      tree->SetBranchAddress("dy1", &directionY) < 0 ||
      tree->SetBranchAddress("dz1", &directionZ) < 0) {
    std::cerr << "Failed to bind scalar branches in tree 'results'." << std::endl;
    tree->ResetBranchAddresses();
    return;
  }
  std::array<std::vector<double>, 14> scalarValues;

  std::vector<double>* clusterX = nullptr;
  std::vector<double>* clusterY = nullptr;
  std::vector<double>* clusterZ = nullptr;
  std::vector<double>* clusterT = nullptr;
  std::vector<double>* clusterE = nullptr;
  std::vector<int>* clusterElectrons = nullptr;
  if (tree->SetBranchAddress("clusterX", &clusterX) < 0 ||
      tree->SetBranchAddress("clusterY", &clusterY) < 0 ||
      tree->SetBranchAddress("clusterZ", &clusterZ) < 0 ||
      tree->SetBranchAddress("clusterT", &clusterT) < 0 ||
      tree->SetBranchAddress("clusterE", &clusterE) < 0 ||
      tree->SetBranchAddress("clusterElectrons", &clusterElectrons) < 0) {
    std::cerr << "Failed to bind cluster branches in tree 'results'." << std::endl;
    tree->ResetBranchAddresses();
    return;
  }

  std::vector<double> clusterDistances3D;
  std::vector<double> clusterTimeDifferences;
  std::vector<double> clusterEnergies;
  std::vector<double> clusterElectronCounts;
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
    if (tree->GetEntry(entry) <= 0) {
      std::cerr << "Failed to read entry " << entry << " from tree 'results'."
                << std::endl;
      tree->ResetBranchAddresses();
      return;
    }
    if (!clusterX || !clusterY || !clusterZ || !clusterT || !clusterE ||
        !clusterElectrons) {
      std::cerr << "A cluster branch is empty at entry " << entry << "."
                << std::endl;
      tree->ResetBranchAddresses();
      return;
    }

    const double amplitude = std::abs(peakAmplitude);
    scalarValues[0].push_back(totalCharge);
    scalarValues[1].push_back(amplitude);
    scalarValues[2].push_back(peakTime);
    scalarValues[3].push_back(riseTime);
    scalarValues[4].push_back(fallTime);
    scalarValues[5].push_back(trackX);
    scalarValues[6].push_back(trackY);
    scalarValues[7].push_back(trackZ);
    scalarValues[8].push_back(directionX);
    scalarValues[9].push_back(directionY);
    scalarValues[10].push_back(directionZ);
    scalarValues[11].push_back(
        totalCharge != 0.0 ? amplitude / totalCharge : 0.0);
    scalarValues[12].push_back(peakTime - riseTime);
    scalarValues[13].push_back(fallTime - riseTime);

    const auto clusterCount = clusterX->size();
    if (clusterY->size() != clusterCount || clusterZ->size() != clusterCount ||
        clusterT->size() != clusterCount || clusterE->size() != clusterCount ||
        clusterElectrons->size() != clusterCount) {
      std::cerr << "Cluster branch sizes do not match at entry " << entry << "."
                << std::endl;
      tree->ResetBranchAddresses();
      return;
    }

    for (std::size_t cluster = 0; cluster < clusterCount; ++cluster) {
      clusterEnergies.push_back(clusterE->at(cluster));
      clusterElectronCounts.push_back(
          static_cast<double>(clusterElectrons->at(cluster)));
      if (cluster + 1 == clusterCount) {
        continue;
      }

      const double dx = clusterX->at(cluster + 1) - clusterX->at(cluster);
      const double dy = clusterY->at(cluster + 1) - clusterY->at(cluster);
      const double dz = clusterZ->at(cluster + 1) - clusterZ->at(cluster);
      clusterDistances3D.push_back(std::sqrt(dx * dx + dy * dy + dz * dz) * 10);
      clusterTimeDifferences.push_back(
          std::abs(clusterT->at(cluster + 1) - clusterT->at(cluster)) * 1000);
    }
  }
  tree->ResetBranchAddresses();

  const auto histogramRange = [](const std::vector<double>& values) {
    if (values.empty()) {
      return std::make_pair(0.0, 1.0);
    }
    const auto bounds = std::minmax_element(values.begin(), values.end());
    double minimum = *bounds.first;
    double maximum = *bounds.second;
    if (minimum == maximum) {
      const double padding = std::max(0.5, std::abs(minimum) * 0.05);
      minimum -= padding;
      maximum += padding;
    } else {
      const double padding = (maximum - minimum) * 0.05;
      minimum -= padding;
      maximum += padding;
    }
    return std::make_pair(minimum, maximum);
  };

  const auto makeHistogram = [&histogramRange](
                                 const char* name, const char* title,
                                 const std::vector<double>& values) {
    const auto range = histogramRange(values);
    auto* histogram = new TH1D(name, title, 100, range.first, range.second);
    for (const double value : values) {
      histogram->Fill(value);
    }
    return histogram;
  };

  outputFile.cd();
  for (std::size_t i = 0; i < specs.size(); ++i) {
    auto* histogram =
        makeHistogram(specs[i].name, specs[i].title, scalarValues[i]);
    if (histogram->GetEntries() != static_cast<double>(scalarValues[i].size())) {
      std::cerr << "Failed to fill all entries in histogram '" << specs[i].name
                << "'." << std::endl;
      return;
    }
    canvas.cd(static_cast<int>(i + 1));
    histogram->Draw("HIST");
    histogram->Write();
  }

  const std::array<TH1D*, 4> clusterHistograms{{
      makeHistogram("h_cluster_distance_3d",
                    "Distance between consecutive clusters;Distance [mm];Pairs",
                    clusterDistances3D),
      makeHistogram("h_cluster_time_difference",
                    "Time difference between consecutive clusters;Time difference [ps];Pairs",
                    clusterTimeDifferences),
      makeHistogram("h_cluster_electrons",
                    "Electrons per cluster;Number of electrons;Clusters",
                    clusterElectronCounts),
      makeHistogram("h_cluster_energy", "Cluster energy;Energy [eV];Clusters",
                    clusterEnergies)}};
  for (std::size_t i = 0; i < clusterHistograms.size(); ++i) {
    auto* pad = canvas.GetPad(static_cast<int>(specs.size() + i + 1));
    pad->cd();
    if (i >= 2) {
      pad->SetLogy();
      clusterHistograms[i]->SetMinimum(0.5);
    }
    clusterHistograms[i]->Draw("HIST");
    clusterHistograms[i]->Write();
  }

  canvas.Write();
  canvas.SaveAs(imagePath.c_str());
  std::cout << "Histograms saved to " << outputPath << std::endl;
  std::cout << "Canvas saved to " << imagePath << std::endl;

  analyse_correlations(runNumber);
}