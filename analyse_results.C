#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <array>
#include <iostream>
#include <string>

void analyse_results(const char* inputPath = "results/results.root",
                     const char* outputPath = "results/results_histograms.root",
                     const char* imagePath = "results/results_histograms.png") {
  TFile inputFile(inputPath, "READ");
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
  const std::array<HistogramSpec, 4> specs{{
      {"totalCharge", "h_totalCharge", "Total charge;Charge [fC];Events"},
      {"q_peak", "h_q_peak", "Peak charge;Peak charge [fC];Events"},
      {"t_peak", "h_t_peak", "Peak time;Time [ns];Events"},
      {"t_threshold", "h_t_threshold", "Threshold crossing time;Time [ns];Events"}}};

  for (const auto& spec : specs) {
    if (!tree->GetBranch(spec.branch)) {
      std::cerr << "Branch '" << spec.branch << "' not found in tree 'results'." << std::endl;
      return;
    }
  }

  TFile outputFile(outputPath, "RECREATE");
  if (outputFile.IsZombie()) {
    std::cerr << "Cannot create output file: " << outputPath << std::endl;
    return;
  }

  outputFile.cd();
  TCanvas canvas("c_results", "Simulation results", 1200, 800);
  canvas.Divide(2, 2);

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
  canvas.SaveAs(imagePath);
  std::cout << "Histograms saved to " << outputPath << std::endl;
  std::cout << "Canvas saved to " << imagePath << std::endl;
}