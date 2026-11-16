#include <TCanvas.h>
#include <TAxis.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Garfield/MediumMagboltz.hh"

R__LOAD_LIBRARY(libGarfield)

using namespace Garfield;

namespace {

struct GasMixture {
  const char* fileName;
  const char* label;
  int color;
};

struct GasCurve {
  std::unique_ptr<TGraph> graph;
  std::string label;
};

bool loadAndSampleGas(const GasMixture& mixture, const std::string& gasDirectory,
                      GasCurve& velocity, GasCurve& longitudinalDiffusion,
                      GasCurve& transverseDiffusion) {
  std::string gasPath = gasDirectory;
  if (!gasPath.empty() && gasPath.back() != '/' && gasPath.back() != '\\') {
    gasPath += '/';
  }
  gasPath += mixture.fileName;

  MediumMagboltz gas;
  if (!gas.LoadGasFile(gasPath)) {
    std::cerr << "Cannot load Garfield gas table: " << gasPath << std::endl;
    return false;
  }

  velocity.graph = std::make_unique<TGraph>();
  velocity.graph->SetName((std::string("velocity_") + mixture.fileName).c_str());
  velocity.label = mixture.label;
  velocity.graph->SetLineColor(mixture.color);
  velocity.graph->SetLineWidth(2);

  longitudinalDiffusion.graph = std::make_unique<TGraph>();
  longitudinalDiffusion.graph->SetName(
      (std::string("diffusionL_") + mixture.fileName).c_str());
  longitudinalDiffusion.label = std::string(mixture.label) + " longitudinale";
  longitudinalDiffusion.graph->SetLineColor(mixture.color);
  longitudinalDiffusion.graph->SetLineStyle(1);
  longitudinalDiffusion.graph->SetLineWidth(2);

  transverseDiffusion.graph = std::make_unique<TGraph>();
  transverseDiffusion.graph->SetName(
      (std::string("diffusionT_") + mixture.fileName).c_str());
  transverseDiffusion.label = std::string(mixture.label) + " trasversale";
  transverseDiffusion.graph->SetLineColor(mixture.color);
  transverseDiffusion.graph->SetLineStyle(2);
  transverseDiffusion.graph->SetLineWidth(2);

  constexpr int nFieldPoints = 200;
  constexpr double eMin = 1.;
  constexpr double eMax = 1.e6;
  const double logEMin = std::log(eMin);
  const double logEMax = std::log(eMax);
  for (int i = 0; i < nFieldPoints; ++i) {
    const double fraction =
        static_cast<double>(i) / static_cast<double>(nFieldPoints - 1);
    const double electricField = std::exp(logEMin + fraction * (logEMax - logEMin));

    double vx = 0.;
    double vy = 0.;
    double vz = 0.;
    if (gas.ElectronVelocity(electricField, 0., 0., 0., 0., 0., vx, vy, vz)) {
      const double speed = std::sqrt(vx * vx + vy * vy + vz * vz);
      if (std::isfinite(speed)) {
        velocity.graph->SetPoint(velocity.graph->GetN(), electricField, speed);
      }
    }

    double dl = 0.;
    double dt = 0.;
    if (gas.ElectronDiffusion(electricField, 0., 0., 0., 0., 0., dl, dt)) {
      if (std::isfinite(dl)) {
        longitudinalDiffusion.graph->SetPoint(
            longitudinalDiffusion.graph->GetN(), electricField, dl);
      }
      if (std::isfinite(dt)) {
        transverseDiffusion.graph->SetPoint(
            transverseDiffusion.graph->GetN(), electricField, dt);
      }
    }
  }

  if (velocity.graph->GetN() < 2 ||
      longitudinalDiffusion.graph->GetN() < 2 ||
      transverseDiffusion.graph->GetN() < 2) {
    std::cerr << "Insufficient Garfield transport data in " << gasPath
              << " over the electric-field range " << eMin << " to " << eMax
              << " V/cm." << std::endl;
    return false;
  }
  return true;
}

void drawComparison(TCanvas& canvas, const std::vector<GasCurve>& curves,
                    TLegend& legend, const char* title,
                    const char* yAxisTitle) {
  double xMinimum = std::numeric_limits<double>::infinity();
  double xMaximum = -std::numeric_limits<double>::infinity();
  double yMaximum = 0.;
  for (const auto& curve : curves) {
    for (int i = 0; i < curve.graph->GetN(); ++i) {
      double x = 0.;
      double y = 0.;
      curve.graph->GetPoint(i, x, y);
      if (std::isfinite(x) && x > 0.) {
        xMinimum = std::min(xMinimum, x);
        xMaximum = std::max(xMaximum, x);
      }
      if (std::isfinite(y)) {
        yMaximum = std::max(yMaximum, y);
      }
    }
  }

  auto* first = curves.front().graph.get();
  first->SetTitle(title);
  first->SetMinimum(0.);
  first->SetMaximum(yMaximum > 0. ? yMaximum * 1.1 : 1.);
  first->Draw("AL");
  first->GetXaxis()->SetTitle("Electric field [V/cm]");
  first->GetYaxis()->SetTitle(yAxisTitle);
  if (std::isfinite(xMinimum) && xMaximum > xMinimum) {
    first->GetXaxis()->SetLimits(xMinimum, xMaximum);
  }
  canvas.SetLogx();
  canvas.SetGrid();

  for (std::size_t i = 1; i < curves.size(); ++i) {
    curves[i].graph->Draw("L SAME");
  }
  legend.SetX1NDC(0.64);
  legend.SetY1NDC(0.68);
  legend.SetX2NDC(0.89);
  legend.SetY2NDC(0.89);
  for (const auto& curve : curves) {
    legend.AddEntry(curve.graph.get(), curve.label.c_str(), "l");
  }
  legend.Draw();
  canvas.Modified();
  canvas.Update();
}

}  // namespace

void analyseGas(const char* gasDirectory = "gas_files",
                const char* outputPath = "results/gasComparison.root") {
  const std::vector<GasMixture> mixtures = {
      {"ar_co2_90_10.gas", "Ar-CO_{2} 90-10", kAzure + 2},
      {"ar_co2_85_15.gas", "Ar-CO_{2} 85-15", kOrange + 7},
      {"ar_co2_70_30.gas", "Ar-CO_{2} 70-30", kViolet + 1}};

  std::vector<GasCurve> velocityCurves;
  std::vector<GasCurve> diffusionCurves;
  for (const auto& mixture : mixtures) {
    GasCurve velocity;
    GasCurve longitudinalDiffusion;
    GasCurve transverseDiffusion;
    if (!loadAndSampleGas(mixture, gasDirectory, velocity,
                          longitudinalDiffusion, transverseDiffusion)) {
      return;
    }
    velocityCurves.push_back(std::move(velocity));
    diffusionCurves.push_back(std::move(longitudinalDiffusion));
    diffusionCurves.push_back(std::move(transverseDiffusion));
  }

  TLegend velocityLegend;
  TLegend diffusionLegend;
  TCanvas velocityComparison("cVelocity", "Electron Drift Velocity", 1000, 700);
  drawComparison(velocityComparison, velocityCurves, velocityLegend,
                 "Electron drift velocity;Electric field [V/cm];"
                 "Drift velocity [cm/ns]",
                 "Drift velocity [cm/ns]");
  TCanvas diffusionComparison("cDiffusion", "Electron Diffusion", 1000, 700);
  drawComparison(diffusionComparison, diffusionCurves, diffusionLegend,
                 "Electron diffusion;Electric field [V/cm];"
                 "Diffusion coefficient [cm^{1/2}]",
                 "Diffusion coefficient [cm^{1/2}]");

  std::unique_ptr<TFile> output(TFile::Open(outputPath, "RECREATE"));
  if (!output || output->IsZombie()) {
    std::cerr << "Cannot create gas-comparison file: " << outputPath
              << std::endl;
    return;
  }
  velocityComparison.Write();
  diffusionComparison.Write();
  velocityComparison.SaveAs("results/gasComparison_velocity.png");
  diffusionComparison.SaveAs("results/gasComparison_diffusion.png");
  output->Close();
  std::cout << "Saved Garfield gas comparison plots to " << outputPath
            << std::endl;
}
