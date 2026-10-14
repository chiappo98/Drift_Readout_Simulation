#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <array>
#include <iostream>
#include <string>

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
    const char* branch;
    const char* name;
    const char* title;
  };
  const std::array<HistogramSpec, 11> specs{{
      {"totalCharge", "h_totalCharge", "Total charge;Charge [fC];Events"},
      {"q_peak", "h_q_peak", "Peak charge;Peak charge [fC];Events"},
      {"t_peak", "h_t_peak", "Peak time;Time [ns];Events"},
      {"t_rise", "h_t_rise", "Rise threshold crossing;Time [ns];Events"},
      {"t_fall", "h_t_fall", "Fall threshold crossing;Time [ns];Events"},
      {"x1", "h_x1", "Track start x;Position [cm];Events"},
      {"y1", "h_y1", "Track start y;Position [cm];Events"},
      {"z1", "h_z1", "Track start z;Position [cm];Events"},
      {"dx1", "h_dx1", "Track direction dx;Direction cosine;Events"},
      {"dy1", "h_dy1", "Track direction dy;Direction cosine;Events"},
      {"dz1", "h_dz1", "Track direction dz;Direction cosine;Events"}}};

  for (const auto& spec : specs) {
    if (!tree->GetBranch(spec.branch)) {
      std::cerr << "Branch '" << spec.branch << "' not found in tree 'results'." << std::endl;
      return;
    }
  }

  TFile outputFile(outputPath.c_str(), "RECREATE");
  if (outputFile.IsZombie()) {
    std::cerr << "Cannot create output file: " << outputPath << std::endl;
    return;
  }

  outputFile.cd();
  TCanvas canvas("c_results", "Simulation results", 1200, 1200);
  canvas.Divide(3, 4);

  for (std::size_t i = 0; i < specs.size(); ++i) {
    const auto& spec = specs[i];
    const std::string expression = std::string(spec.branch) + ">>" + spec.name + "(100)";
    if (tree->Draw(expression.c_str(), "", "goff") < 0) {
      std::cerr << "Failed to fill histogram for branch '" << spec.branch << "'." << std::endl;
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

  canvas.Write();
  canvas.SaveAs(imagePath.c_str());
  std::cout << "Histograms saved to " << outputPath << std::endl;
  std::cout << "Canvas saved to " << imagePath << std::endl;
}