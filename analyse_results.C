#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TPad.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

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
      {"abs(q_peak)", "h_q_peak", "Peak amplitude;Amplitude [shaper units];Events"},
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
       "Fall time - rise time;Time difference [ns];Events"}}};

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

  for (std::size_t i = 0; i < specs.size(); ++i) {
    const auto& spec = specs[i];
    const std::string expression = std::string(spec.expression) + ">>" + spec.name + "(100)";
    if (tree->Draw(expression.c_str(), "", "goff") < 0) {
      std::cerr << "Failed to fill histogram '" << spec.name << "'." << std::endl;
      return;
    }

    auto* histogram = dynamic_cast<TH1*>(outputFile.Get(spec.name));
    if (!histogram) {
      std::cerr << "Failed to create histogram '" << spec.name << "'." << std::endl;
      return;
    }
    histogram->SetTitle(spec.title);
    canvas.cd(static_cast<int>(i + 1));
    histogram->Draw("HIST");
    histogram->Write();
  }

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
}