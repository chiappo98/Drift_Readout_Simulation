#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TLegend.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace {

double transfer1Order(double* x, double* parameters) {
  const double scaledTime = x[0] / parameters[0];
  return scaledTime * std::exp(1.0 - scaledTime);
}

double transfer3Order(double* x, double* parameters) {
  const double scaledTime = x[0] / parameters[0];
  return std::pow(scaledTime, 3) * std::exp(3.0 * (1.0 - scaledTime));
}

double transferVtx(double* x, double*) {
  constexpr double gain = 39.6 / 40.0;
  constexpr double order = 3.34;
  constexpr double tauNs = 3.06;
  constexpr double offsetNs = 1.463;
  const double timeSinceOffset = x[0] - offsetNs;
  if (timeSinceOffset <= 0.0) {
    return 0.0;
  }

  const double response =
      gain * std::pow(timeSinceOffset / (order * tauNs), order) *
      std::exp(order - timeSinceOffset / tauNs);
  return response / std::abs(gain);
}

}  // namespace

void plot_transfer_functions(double shapingTimeNs = 25.0,
                             const char* outputPath =
                                 "plot/transfer_functions.png") {
  if (!std::isfinite(shapingTimeNs) || shapingTimeNs <= 0.0) {
    std::cerr << "The shaping time must be a positive finite value." << std::endl;
    return;
  }
  if (!outputPath || outputPath[0] == '\0') {
    std::cerr << "An output path is required." << std::endl;
    return;
  }

  constexpr double vtxOffsetNs = 1.463;
  constexpr double vtxTauNs = 3.06;
  constexpr int functionSamples = 300;
  const double endTimeNs =
      std::max(5.0 * shapingTimeNs, vtxOffsetNs + 8.0 * vtxTauNs);

  auto* order1 = new TF1("transfer_order_1", transfer1Order, 0.0, endTimeNs, 1);
  order1->SetParameter(0, shapingTimeNs);
  order1->SetNpx(functionSamples);
  order1->SetTitle(
      "VTX and shaper transfer functions;Time [ns];Normalized response");
  order1->SetLineColor(kBlue + 1);
  order1->SetLineWidth(2);

  auto* order3 = new TF1("transfer_order_3", transfer3Order, 0.0, endTimeNs, 1);
  order3->SetParameter(0, shapingTimeNs);
  order3->SetNpx(functionSamples);
  order3->SetLineColor(kRed + 1);
  order3->SetLineWidth(2);

  auto* vtx = new TF1("transfer_vtx", transferVtx, 0.0, endTimeNs, 0);
  vtx->SetNpx(functionSamples);
  vtx->SetLineColor(kGreen + 2);
  vtx->SetLineWidth(2);

  auto* canvas = new TCanvas("c_transfer_functions",
                             "VTX and shaper transfer functions", 1000, 700);
  order1->SetMinimum(-0.15);
  order1->SetMaximum(1.1);
  order1->Draw();
  order3->Draw("SAME");
  vtx->Draw("SAME");

  auto* legend = new TLegend(0.58, 0.68, 0.88, 0.88);
  legend->AddEntry(vtx, "VTX", "l");
  const std::string order1Label =
      "1st order - #tau = " + std::to_string((int)shapingTimeNs) + " ns";
  legend->AddEntry(order1, order1Label.c_str(), "l");
  const std::string order3Label =
      "3rd order - #tau = " + std::to_string((int)shapingTimeNs) + " ns";
  legend->AddEntry(order3, order3Label.c_str(), "l");
  legend->Draw();

  canvas->SaveAs(outputPath);
  std::ifstream outputFile(outputPath, std::ios::binary);
  if (!outputFile) {
    std::cerr << "Failed to create transfer-function plot: " << outputPath
              << std::endl;
    return;
  }
  std::cout << "Transfer functions saved to " << outputPath << std::endl;
}
