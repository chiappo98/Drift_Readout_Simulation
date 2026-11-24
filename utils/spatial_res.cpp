#include <TAxis.h>
#include <array>
#include <TCanvas.h>
#include <TCollection.h>
#include <TColor.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TSpline.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int kSlices = 23;
constexpr double kSplineMinimumTrackY = -0.65;
constexpr double kSplineMaximumTrackY = 0.0;
constexpr double kResolutionFitMinimumTrackY = -0.6;
constexpr int kResidualBins = 14;

double crystalBall(double* x, double* parameters) {
  const double alpha = std::abs(parameters[3]);
  const double exponent = parameters[4];
  if (alpha <= 0.0 || exponent <= 0.0 || parameters[2] <= 0.0) {
    return 0.0;
  }

  const double standardized = (x[0] - parameters[1]) / parameters[2];
  if (standardized < alpha) {
    return parameters[0] * std::exp(-0.5 * standardized * standardized);
  }

  const double coefficient =
      std::pow(exponent / alpha, exponent) * std::exp(-0.5 * alpha * alpha);
  const double offset = exponent / alpha - alpha;
  return parameters[0] * coefficient *
         std::pow(offset + standardized, -exponent);
}

double crystalBallWithLowSideOffset(double* x, double* parameters) {
  const double lowSideOffset = x[0] < parameters[1] ? parameters[5] : 0.0;
  return crystalBall(x, parameters) + lowSideOffset;
}

double crystalBallLeftTail(double* x, double* parameters) {
  double reflectedX = 2.0 * parameters[1] - x[0];
  return crystalBall(&reflectedX, parameters);
}

double residualResolutionComponents(double* x, double* parameters) {
  const double primaryIonization =
      parameters[0] * std::exp(-parameters[0] * (-x[0])); //std::abs
  const double electronics =
      parameters[2] + parameters[3] * std::exp(parameters[4] * x[0]);
  if (!std::isfinite(primaryIonization) || !std::isfinite(electronics)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double variance = primaryIonization * primaryIonization +
                          parameters[1] * (-x[0]) +    //std::abs
                          electronics * electronics;
  return variance > 0.0 && std::isfinite(variance)
             ? std::sqrt(variance)
             : std::numeric_limits<double>::quiet_NaN();
}

double primaryIonizationResolution(double* x, double* parameters) {
  return parameters[0] * std::exp(-parameters[0] *(-x[0]));  // std::abs
}

double diffusionResolution(double* x, double* parameters) {
  return std::sqrt(std::max(0.0, parameters[0] * (-x[0]))); //std::abs(x[0])
}

double electronicsResolution(double* x, double* parameters) {
  const double component =
      parameters[0] + parameters[1] * std::exp(parameters[2] * x[0]);
  return std::isfinite(component) ? component
                                  : std::numeric_limits<double>::quiet_NaN();
}

std::string twoDigits(int value) {
  std::ostringstream result;
  result << std::setw(2) << std::setfill('0') << value;
  return result.str();
}

double readConfigNumber(const std::string& configSnapshotPath,
                        const std::string& key) {
  std::ifstream input(configSnapshotPath);
  if (!input) {
    throw std::runtime_error("Cannot open run configuration snapshot: " +
                             configSnapshotPath);
  }
  const std::string contents((std::istreambuf_iterator<char>(input)),
                             std::istreambuf_iterator<char>());
  const std::regex pattern("\"" + key +
                           "\"\\s*:\\s*(-?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)"
                           "(?:[eE][+-]?[0-9]+)?)");
  std::smatch match;
  if (!std::regex_search(contents, match, pattern)) {
    throw std::runtime_error("Missing numeric '" + key +
                             "' in configuration snapshot " +
                             configSnapshotPath);
  }
  return std::stod(match[1].str());
}

std::string formatFieldFilenameComponent(double value, int precision = 3) {
  std::ostringstream stream;
  stream << std::fixed << std::setprecision(precision) << std::abs(value);
  std::string result = stream.str();
  while (result.size() > 1 && result.back() == '0') {
    result.pop_back();
  }
  if (!result.empty() && result.back() == '.') {
    result.pop_back();
  }
  return result;
}

std::string fieldPlotPathFromSnapshot(const std::string& runBase) {
  const std::string snapshotPath = runBase + "_configurations.txt";
  const double argon = readConfigNumber(snapshotPath, "argon_percent");
  const double co2 = readConfigNumber(snapshotPath, "co2_percent");
  const double field = readConfigNumber(snapshotPath, "vField");
  const double mylar = readConfigNumber(snapshotPath, "vMylar");
  const double sense = readConfigNumber(snapshotPath, "vSense");
  const double cellW = readConfigNumber(snapshotPath, "cell_w");
  const double cellY = readConfigNumber(snapshotPath, "cell_y");
  const double dSense = readConfigNumber(snapshotPath, "dSense");
  const double dField = readConfigNumber(snapshotPath, "dField");
  return "results/Ar-" + formatFieldFilenameComponent(argon, 0) +
         "_CO2-" + formatFieldFilenameComponent(co2, 0) +
         "_vField-" + formatFieldFilenameComponent(field) +
         "_vMylar-" + formatFieldFilenameComponent(mylar) +
         "_vSense-" + formatFieldFilenameComponent(sense) +
         "_cellw-" + formatFieldFilenameComponent(cellW) +
         "_celly-" + formatFieldFilenameComponent(cellY) +
         "_dSense-" + formatFieldFilenameComponent(dSense) +
         "_dField-" + formatFieldFilenameComponent(dField) + ".root";
}

void collectGraphs(TVirtualPad& pad, std::vector<TGraph*>& graphs) {
  TIter next(pad.GetListOfPrimitives());
  while (auto* object = next()) {
    if (auto* graph = dynamic_cast<TGraph*>(object)) {
      if (graph->GetN() >= 2) {
        graphs.push_back(graph);
      }
    } else if (auto* multiGraph = dynamic_cast<TMultiGraph*>(object)) {
      TIter nextGraph(multiGraph->GetListOfGraphs());
      while (auto* graphObject = nextGraph()) {
        auto* graph = dynamic_cast<TGraph*>(graphObject);
        if (graph && graph->GetN() >= 2) {
          graphs.push_back(graph);
        }
      }
    } else if (auto* nestedPad = dynamic_cast<TVirtualPad*>(object)) {
      collectGraphs(*nestedPad, graphs);
    }
  }
}

bool interpolateGraph(const TGraph& graph, double x, double& y) {
  std::vector<std::pair<double, double>> points;
  points.reserve(graph.GetN());
  for (int point = 0; point < graph.GetN(); ++point) {
    double graphX = 0.0;
    double graphY = 0.0;
    graph.GetPoint(point, graphX, graphY);
    if (std::isfinite(graphX) && std::isfinite(graphY)) {
      points.emplace_back(graphX, graphY);
    }
  }
  std::sort(points.begin(), points.end());
  if (points.size() < 2 || x < points.front().first ||
      x > points.back().first) {
    return false;
  }

  const auto upper = std::lower_bound(
      points.begin(), points.end(), x,
      [](const std::pair<double, double>& point, double value) {
        return point.first < value;
      });
  if (upper == points.begin()) {
    y = upper->second;
    return true;
  }
  if (upper == points.end()) {
    y = points.back().second;
    return true;
  }
  if (upper->first == x) {
    y = upper->second;
    return true;
  }
  const auto lower = upper - 1;
  const double span = upper->first - lower->first;
  if (span <= 0.0) {
    return false;
  }
  const double fraction = (x - lower->first) / span;
  y = lower->second + fraction * (upper->second - lower->second);
  return std::isfinite(y);
}

int makeGarfieldVelocityCurves(const std::string& fieldPlotPath,
                               std::array<TGraph, 3>& velocityCurves) {
  const std::string gasPlotPath = "results/gasParameters.root";
  std::unique_ptr<TFile> gasFile(TFile::Open(gasPlotPath.c_str(), "READ"));
  if (!gasFile || gasFile->IsZombie()) {
    std::cerr << "Cannot open Garfield gas-parameter file: " << gasPlotPath
              << std::endl;
    return 0;
  }
  std::unique_ptr<TFile> fieldFile(TFile::Open(fieldPlotPath.c_str(), "READ"));
  if (!fieldFile || fieldFile->IsZombie()) {
    std::cerr << "Cannot open Garfield field-plot file: " << fieldPlotPath
              << std::endl;
    return 0;
  }

  auto* gasCanvas = dynamic_cast<TCanvas*>(gasFile->Get("cm"));
  auto* fieldCanvas = dynamic_cast<TCanvas*>(fieldFile->Get("cf1"));
  if (!gasCanvas || !fieldCanvas) {
    std::cerr << "Required canvas 'cm' or 'cf1' is missing from Garfield "
                 "plot files."
              << std::endl;
    return 0;
  }
  auto* fieldHistogram =
      dynamic_cast<TH2*>(fieldCanvas->GetPrimitive("hE"));
  if (!fieldHistogram) {
    std::cerr << "Electric-field histogram 'hE' is missing from canvas 'cf1'."
              << std::endl;
    return 0;
  }

  constexpr std::array<const char*, 3> labels = {
      "Ar-CO_{2} 90-10", "Ar-CO_{2} 85-15", "Ar-CO_{2} 70-30"};
  constexpr std::array<int, 3> colors = {kCyan + 2, kOrange + 7, kViolet + 1};
  const double xMinimum = fieldHistogram->GetXaxis()->GetXmin();
  const double xMaximum = fieldHistogram->GetXaxis()->GetXmax();
  if (xMinimum > 0.0 || xMaximum < 0.0) {
    std::cerr << "Field histogram does not cover x = 0." << std::endl;
    return 0;
  }
  const int fieldXBin = fieldHistogram->GetXaxis()->FindFixBin(0.0);
  int fieldXLowBin = fieldXBin;
  int fieldXHighBin = fieldXBin;
  if (fieldHistogram->GetXaxis()->GetBinCenter(fieldXBin) < 0.0) {
    fieldXHighBin = fieldXBin + 1;
  } else if (fieldHistogram->GetXaxis()->GetBinCenter(fieldXBin) > 0.0) {
    fieldXLowBin = fieldXBin - 1;
  }
  const double xLow = fieldHistogram->GetXaxis()->GetBinCenter(fieldXLowBin);
  const double xHigh = fieldHistogram->GetXaxis()->GetBinCenter(fieldXHighBin);
  if (xLow > 0.0 || xHigh < 0.0) {
    std::cerr << "Cannot interpolate the field histogram at x = 0."
              << std::endl;
    return 0;
  }

  int populatedCurves = 0;
  for (int gasIndex = 0; gasIndex < 3; ++gasIndex) {
    auto* gasPad = gasCanvas->GetPad(gasIndex + 1);
    if (!gasPad) {
      std::cerr << "Gas-parameter canvas 'cm' is missing pad "
                << gasIndex + 1 << "." << std::endl;
      continue;
    }
    std::vector<TGraph*> graphs;
    collectGraphs(*gasPad, graphs);
    const auto graph = std::max_element(
        graphs.begin(), graphs.end(),
        [](const TGraph* left, const TGraph* right) {
          return left->GetN() < right->GetN();
        });
    if (graph == graphs.end()) {
      std::cerr << "No velocity graph found in gas-parameter pad "
                << gasIndex + 1 << "." << std::endl;
      continue;
    }

    TGraph& velocityCurve = velocityCurves[gasIndex];
    velocityCurve.SetName(
        ("garfield_drift_velocity_ArCO2_" +
         std::to_string(gasIndex + 1))
            .c_str());
    velocityCurve.SetTitle(
        "Garfield drift velocity vs track y;Track y [cm];"
        "Drift velocity [cm/#mus]");
    velocityCurve.SetLineColor(colors[gasIndex]);
    velocityCurve.SetLineWidth(2);
    for (int yBin = 1; yBin <= fieldHistogram->GetNbinsY(); ++yBin) {
      const double trackY = fieldHistogram->GetYaxis()->GetBinCenter(yBin);
      if (trackY < -0.7 || trackY > 0.0) {
        continue;
      }
      const double lowField =
          fieldHistogram->GetBinContent(fieldXLowBin, yBin);
      const double highField =
          fieldHistogram->GetBinContent(fieldXHighBin, yBin);
      if (!std::isfinite(lowField) || !std::isfinite(highField)) {
        continue;
      }
      double electricField = lowField;
      if (fieldXHighBin != fieldXLowBin && xHigh != xLow) {
        const double fraction = (0.0 - xLow) / (xHigh - xLow);
        electricField = lowField + fraction * (highField - lowField);
      }
      double velocity = 0.0;
      if (std::isfinite(electricField) && electricField >= 0.0 &&
          interpolateGraph(**graph, electricField, velocity)) {
        velocityCurve.SetPoint(velocityCurve.GetN(), trackY, velocity * 1000.0);
      }
    }
    if (velocityCurve.GetN() > 0) {
      ++populatedCurves;
    } else {
      std::cerr << "No Garfield velocity samples overlap the electric-field "
                   "range for "
                << labels[gasIndex] << "." << std::endl;
    }
  }
  return populatedCurves;
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
                         crystalBallWithLowSideOffset,
                         distribution->GetXaxis()->GetXmin(),
                         distribution->GetXaxis()->GetXmax(), 6);
      crystalBallFit.SetParNames("Amplitude", "MPV", "Sigma", "Alpha", "N",
                                 "Low-side offset");
      crystalBallFit.SetParameters(distribution->GetMaximum(), meanStart,
                                   sigmaStart, 1.5, 3.0,
                                   distribution->GetMinimum());
      crystalBallFit.SetParLimits(0, 0.0, distribution->GetMaximum() * 100.0);
      crystalBallFit.SetParLimits(1, distribution->GetXaxis()->GetXmin(),
                                  distribution->GetXaxis()->GetXmax());
      crystalBallFit.SetParLimits(2, binWidth * 0.1,
                                  distribution->GetXaxis()->GetXmax() -
                                      distribution->GetXaxis()->GetXmin());
      crystalBallFit.SetParLimits(3, 0.1, 10.0);
      crystalBallFit.SetParLimits(4, 1.01, 100.0);
      crystalBallFit.SetParLimits(5, 0.0, distribution->GetMaximum());
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
  TGraph splineVelocityCurve;
  std::unique_ptr<TF1> mpvLineFit;
  TGraph lineCurve;
  TGraph lineVelocityCurve;
  splineCurve.SetName("crystal_ball_mpv_spline_curve");
  splineCurve.SetLineColor(kGreen + 2);
  splineCurve.SetLineWidth(3);
  splineVelocityCurve.SetName("spline_drift_velocity_vs_trackY");
  splineVelocityCurve.SetTitle(
      "Drift velocity vs track y;Track y [cm];Drift velocity [cm/#mus]");
  splineVelocityCurve.SetLineColor(kGreen + 2);
  splineVelocityCurve.SetLineWidth(3);
  lineCurve.SetName("crystal_ball_mpv_linear_fit_curve");
  lineCurve.SetLineColor(kMagenta + 1);
  lineCurve.SetLineWidth(3);
  lineVelocityCurve.SetName("linear_drift_velocity_vs_trackY");
  lineVelocityCurve.SetTitle(
      "Drift velocity vs track y;Track y [cm];Drift velocity [cm/#mus]");
  lineVelocityCurve.SetLineColor(kMagenta + 1);
  lineVelocityCurve.SetLineStyle(2);
  lineVelocityCurve.SetLineWidth(3);
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
    constexpr int kSplineSamples = 300;
    mpvSpline = std::make_unique<TSpline3>(
        "crystal_ball_mpv_spline", &splinePoints, "b2e2");
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

    mpvLineFit = std::make_unique<TF1>(
        "crystal_ball_mpv_linear_fit", "pol1", splineMinimum, splineMaximum);
    mpvLineFit->SetParNames("Intercept", "Slope");
    const int fitStatus = mpvGraph.Fit(mpvLineFit.get(), "QRS0");
    const double slope = mpvLineFit->GetParameter(1);
    if (fitStatus != 0 || !std::isfinite(slope) || slope == 0.0) {
      std::cerr << "Linear MPV fit failed or has an invalid slope (status "
                << fitStatus << ")." << std::endl;
      mpvLineFit.reset();
    } else {
      for (int sample = 0; sample < kSplineSamples; ++sample) {
        const double x =
            splineMinimum +
            (splineMaximum - splineMinimum) * sample /
                (kSplineSamples - 1);
        lineCurve.SetPoint(sample, x, mpvLineFit->Eval(x));
      }
      lineCurve.Draw("L SAME");

      const double lineVelocity = 1000.0 / std::abs(slope);
      if (std::isfinite(lineVelocity)) {
        lineVelocityCurve.SetPoint(0, splineMinimum, lineVelocity);
        lineVelocityCurve.SetPoint(1, splineMaximum, lineVelocity);
      }
    }

    for (int sample = 0; sample < kSplineSamples; ++sample) {
      const double x =
          splineMinimum +
          (splineMaximum - splineMinimum) * sample /
              (kSplineSamples - 1);
      const double timeDerivative = mpvSpline->Derivative(x);
      if (std::isfinite(timeDerivative) && timeDerivative != 0.0) {
        const double velocity = 1000.0 / std::abs(timeDerivative);
        if (std::isfinite(velocity)) {
          splineVelocityCurve.SetPoint(splineVelocityCurve.GetN(), x,
                                       velocity);
        }
      }
    }
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
  if (mpvLineFit) {
    legend.AddEntry(&lineCurve, "MPV linear fit", "l");
  }
  legend.Draw();
  if (mpvGraph.Write() <= 0 ||
      (mpvSpline &&
       (mpvSpline->Write() <= 0 || splineCurve.Write() <= 0)) ||
      (mpvLineFit &&
       (mpvLineFit->Write() <= 0 || lineCurve.Write() <= 0)) ||
      canvas.Write() <= 0) {
    std::cerr << "Failed to write spatial-resolution output: " << rootPath
              << std::endl;
    return;
  }
  canvas.SaveAs(plotPath.c_str());

  if (!mpvSpline) {
    std::cerr << "Cannot calculate residuals or drift velocity without a "
                 "valid MPV spline."
              << std::endl;
    return;
  }

  TMultiGraph driftVelocity;
  driftVelocity.SetName("drift_velocity_vs_trackY");
  driftVelocity.SetTitle(
      "Drift velocity vs track y;Track y [cm];Drift velocity [cm/#mus]");
  if (splineVelocityCurve.GetN() > 0) {
    driftVelocity.Add(&splineVelocityCurve, "L");
  } else {
    std::cerr << "Cannot calculate spline drift velocity in the fit range."
              << std::endl;
  }
  if (mpvLineFit) {
    driftVelocity.Add(&lineVelocityCurve, "L");
  }
  std::array<TGraph, 3> garfieldVelocityCurves;
  int garfieldCurveCount = 0;
  try {
    garfieldCurveCount = makeGarfieldVelocityCurves(
        fieldPlotPathFromSnapshot(runBase), garfieldVelocityCurves);
  } catch (const std::exception& error) {
    std::cerr << "Cannot load Garfield drift-velocity curves: " << error.what()
              << std::endl;
  }
  constexpr std::array<const char*, 3> garfieldGasLabels = {
      "Garfield Ar-CO_{2} 90-10", "Garfield Ar-CO_{2} 85-15",
      "Garfield Ar-CO_{2} 70-30"};
  for (int gasIndex = 0; gasIndex < 3; ++gasIndex) {
    if (garfieldVelocityCurves[gasIndex].GetN() > 0) {
      driftVelocity.Add(&garfieldVelocityCurves[gasIndex], "L");
    }
  }
  TCanvas driftVelocityCanvas("c_drift_velocity",
                              "Drift velocity vs track y", 1000, 700);
  driftVelocity.Draw("A");
  TLegend driftVelocityLegend(0.65, 0.78, 0.88, 0.89);
  if (splineVelocityCurve.GetN() > 0) {
    driftVelocityLegend.AddEntry(&splineVelocityCurve,
                                 "Spline-derived velocity", "l");
  }
  if (mpvLineFit) {
    driftVelocityLegend.AddEntry(&lineVelocityCurve, "Linear-fit velocity",
                                 "l");
  }
  for (int gasIndex = 0; gasIndex < 3; ++gasIndex) {
    if (garfieldVelocityCurves[gasIndex].GetN() > 0) {
      driftVelocityLegend.AddEntry(&garfieldVelocityCurves[gasIndex],
                                   garfieldGasLabels[gasIndex], "l");
    }
  }
  driftVelocityLegend.Draw();
  outputFile->cd();
  if ((splineVelocityCurve.GetN() > 0 &&
       splineVelocityCurve.Write() <= 0) ||
      (mpvLineFit && lineVelocityCurve.Write() <= 0) ||
      driftVelocity.Write() <= 0 || driftVelocityCanvas.Write() <= 0) {
    std::cerr << "Failed to write drift-velocity output to " << rootPath
              << std::endl;
    return;
  }
  for (const auto& velocityCurve : garfieldVelocityCurves) {
    if (velocityCurve.GetN() > 0 && velocityCurve.Write() <= 0) {
      std::cerr << "Failed to write a Garfield drift-velocity curve to "
                << rootPath << std::endl;
      return;
    }
  }
  const std::string driftVelocityPlotPath =
      runBase + "_drift_velocity.png";
  driftVelocityCanvas.SaveAs(driftVelocityPlotPath.c_str());

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
  std::array<std::unique_ptr<TH1D>, kResidualBins>
      linearResidualHistograms;
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
    linearResidualHistograms[bin] = std::make_unique<TH1D>(
        ("trackY_linear_residual_bin_" + twoDigits(bin + 1)).c_str(),
        (title.str() + " (linear fit);Reconstructed y - true y [cm];Events")
            .c_str(),
        390, -0.65, 0.65);
    linearResidualHistograms[bin]->SetDirectory(nullptr);
    linearResidualHistograms[bin]->SetLineColor(kMagenta + 1);
    linearResidualHistograms[bin]->SetLineWidth(2);
    linearResidualHistograms[bin]->SetFillStyle(0);
  }

  Long64_t reconstructedEntries = 0;
  Long64_t rejectedEntries = 0;
  Long64_t linearReconstructedEntries = 0;
  Long64_t linearRejectedEntries = 0;
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

    int bin = static_cast<int>(
        (trackY - kSplineMinimumTrackY) / residualBinWidth);
    if (bin == kResidualBins) {
      bin = kResidualBins - 1;
    }

    double reconstructedTrackY = 0.0;
    if (invertSpline(*mpvSpline, splineMinimum, splineMaximum, riseTime,
                     reconstructedTrackY)) {
      residualHistograms[bin]->Fill(reconstructedTrackY - trackY);
      ++reconstructedEntries;
    } else {
      ++rejectedEntries;
    }

    if (mpvLineFit) {
      const double linearReconstructedTrackY =
          (riseTime - mpvLineFit->GetParameter(0)) /
          mpvLineFit->GetParameter(1);
      if (!std::isfinite(linearReconstructedTrackY) ||
          linearReconstructedTrackY < splineMinimum ||
          linearReconstructedTrackY > splineMaximum) {
        ++linearRejectedEntries;
      } else {
        linearResidualHistograms[bin]->Fill(linearReconstructedTrackY -
                                            trackY);
        ++linearReconstructedEntries;
      }
    }
  }
  eventTree->ResetBranchAddresses();

  TCanvas residualCanvas("c_trackY_residuals",
                         "Track y residuals by true position", 1500, 900);
  residualCanvas.Divide(4, 4);
  TGraphErrors residualSigmaGraph;
  residualSigmaGraph.SetName("residual_crystal_ball_sigma_vs_slice");
  residualSigmaGraph.SetTitle(
      "Residual width vs true track y;True track y slice center [cm];"
      "Width [#mum]");
  residualSigmaGraph.SetMarkerStyle(20);
  residualSigmaGraph.SetMarkerColor(kBlue + 1);
  residualSigmaGraph.SetLineColor(kBlue + 1);
  TGraphErrors linearResidualSigmaGraph;
  linearResidualSigmaGraph.SetName(
      "linear_residual_crystal_ball_sigma_vs_slice");
  linearResidualSigmaGraph.SetTitle(
      "Linear-fit residual width vs true track y;"
      "True track y slice center [cm];Width [#mum]");
  linearResidualSigmaGraph.SetMarkerStyle(22);
  linearResidualSigmaGraph.SetMarkerColor(kMagenta + 1);
  linearResidualSigmaGraph.SetLineColor(kMagenta + 1);
  TGraphErrors residualMpvGraph;
  residualMpvGraph.SetName("residual_crystal_ball_mpv_vs_slice");
  residualMpvGraph.SetTitle(
      "Spline residual MPV vs true track y;"
      "True track y slice center [cm];Crystal Ball MPV [#mum]");
  residualMpvGraph.SetMarkerStyle(20);
  residualMpvGraph.SetMarkerColor(kBlue + 1);
  residualMpvGraph.SetLineColor(kBlue + 1);
  TGraphErrors linearResidualMpvGraph;
  linearResidualMpvGraph.SetName("linear_residual_crystal_ball_mpv_vs_slice");
  linearResidualMpvGraph.SetTitle(
      "Linear-fit residual MPV vs true track y;"
      "True track y slice center [cm];Crystal Ball MPV [#mum]");
  linearResidualMpvGraph.SetMarkerStyle(22);
  linearResidualMpvGraph.SetMarkerColor(kMagenta + 1);
  linearResidualMpvGraph.SetLineColor(kMagenta + 1);
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
        std::make_unique<TLegend>(0.47, 0.70, 0.89, 0.89);
    const auto fitResidualHistogram =
        [&](TH1D& histogram, const std::string& fitName, int color,
            const char* label, TGraphErrors& sigmaGraph,
            TGraphErrors& mpvGraph) {
          if (histogram.GetEntries() <= 0.0 || histogram.GetMaximum() <= 0.0) {
            std::cerr << label << " residual slice " << bin + 1
                      << " is empty; skipping Crystal Ball fit." << std::endl;
            return;
          }
          const double binWidth = histogram.GetXaxis()->GetBinWidth(1);
          const double sigmaStart = std::max(histogram.GetRMS(), binWidth);
          const double meanStart = histogram.GetXaxis()->GetBinCenter(
              histogram.GetMaximumBin());
          TF1 residualFit(fitName.c_str(), crystalBallLeftTail, -0.3, 0.3, 5);
          residualFit.SetParNames("Amplitude", "MPV", "Sigma", "Alpha", "N");
          residualFit.SetParameters(histogram.GetMaximum(), meanStart,
                                    sigmaStart, 1.5, 3.0);
          residualFit.SetParLimits(0, 0.0, histogram.GetMaximum() * 100.0);
          residualFit.SetParLimits(1, -0.3, 0.3);
          residualFit.SetParLimits(2, binWidth * 0.1, 1.3);
          residualFit.SetParLimits(3, 0.1, 10.0);
          residualFit.SetParLimits(4, 1.01, 100.0);
          const int fitStatus = histogram.Fit(&residualFit, "QRS0");
          if (fitStatus != 0) {
            std::cerr << "Crystal Ball fit failed for " << label
                      << " residual slice " << bin + 1 << " (status "
                      << fitStatus << ")." << std::endl;
            return;
          }
          auto* fittedFunction = histogram.GetFunction(fitName.c_str());
          if (!fittedFunction) {
            std::cerr << "Crystal Ball function missing for " << label
                      << " residual slice " << bin + 1 << "." << std::endl;
            return;
          }
          fittedFunction->SetLineColor(color);
          fittedFunction->SetLineWidth(3);
          residualLegends[bin]->AddEntry(&histogram, label, "l");
          residualLegends[bin]->AddEntry(fittedFunction,
                                         "Crystal Ball fit", "l");
          sigmaGraph.SetPoint(sigmaGraph.GetN(), sliceCenter,
                              residualFit.GetParameter(2) * 10000.0);
          sigmaGraph.SetPointError(sigmaGraph.GetN() - 1, 0.0,
                                   residualFit.GetParError(2) * 10000.0);
          mpvGraph.SetPoint(mpvGraph.GetN(), sliceCenter,
                            residualFit.GetParameter(1) * 10000.0);
          mpvGraph.SetPointError(mpvGraph.GetN() - 1, 0.0,
                                 residualFit.GetParError(1) * 10000.0);
        };
    fitResidualHistogram(
        residualHistogram,
        "residual_crystal_ball_fit_" + twoDigits(bin + 1), kRed + 1,
        "Spline residuals", residualSigmaGraph, residualMpvGraph);
    fitResidualHistogram(
        *linearResidualHistograms[bin],
        "linear_residual_crystal_ball_fit_" + twoDigits(bin + 1),
        kMagenta + 1, "Linear residuals", linearResidualSigmaGraph,
        linearResidualMpvGraph);

    outputFile->cd();
    if (residualHistogram.Write() <= 0 ||
        linearResidualHistograms[bin]->Write() <= 0) {
      std::cerr << "Failed to write residual histogram for trackY bin "
                << bin + 1 << "." << std::endl;
      return;
    }
    residualCanvas.cd(bin + 1);
    residualHistogram.Draw("HIST");
    linearResidualHistograms[bin]->Draw("HIST SAME");
    if (residualHistogram.GetFunction(
            ("residual_crystal_ball_fit_" + twoDigits(bin + 1)).c_str())) {
      auto* fittedFunction = residualHistogram.GetFunction(
          ("residual_crystal_ball_fit_" + twoDigits(bin + 1)).c_str());
      fittedFunction->Draw("SAME");
    }
    if (linearResidualHistograms[bin]->GetFunction(
            ("linear_residual_crystal_ball_fit_" + twoDigits(bin + 1))
                .c_str())) {
      auto* fittedFunction = linearResidualHistograms[bin]->GetFunction(
          ("linear_residual_crystal_ball_fit_" + twoDigits(bin + 1))
              .c_str());
      fittedFunction->Draw("SAME");
    }
    residualLegends[bin]->Draw();
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

  TMultiGraph residualMeans;
  residualMeans.SetName("residual_mpv_vs_slice");
  residualMeans.SetTitle(
      "Crystal Ball residual MPV vs true track y;"
      "True track y slice center [cm];Crystal Ball MPV [#mum]");
  residualMeans.Add(&residualMpvGraph, "P");
  residualMeans.Add(&linearResidualMpvGraph, "P");
  TCanvas meansCanvas("c_residual_mpv",
                      "Crystal Ball residual MPV by true position", 1000,
                      700);
  residualMeans.Draw("A");
  residualMeans.GetXaxis()->SetLimits(kSplineMinimumTrackY,
                                      kSplineMaximumTrackY);
  residualMeans.GetYaxis()->SetTitle("Crystal Ball MPV [#mum]");
  TLegend meansLegend(0.64, 0.78, 0.88, 0.89);
  meansLegend.AddEntry(&residualMpvGraph, "Spline residuals", "p");
  meansLegend.AddEntry(&linearResidualMpvGraph, "Linear-fit residuals", "p");
  meansLegend.Draw();
  outputFile->cd();
  if (residualMpvGraph.Write() <= 0 ||
      linearResidualMpvGraph.Write() <= 0 || residualMeans.Write() <= 0 ||
      meansCanvas.Write() <= 0) {
    std::cerr << "Failed to write residual-MPV output to " << rootPath
              << std::endl;
    return;
  }

  residualRmsGraph.SetLineWidth(2);
  TMultiGraph residualWidths;
  residualWidths.SetName("residual_widths_vs_slice");
  residualWidths.SetTitle(
      "Residual widths by true track y;True track y slice center [cm];"
      "Width [#mum]");
  residualWidths.Add(&residualSigmaGraph, "P");
  residualWidths.Add(&linearResidualSigmaGraph, "P");
  residualWidths.Add(&residualRmsGraph, "LP");
  TCanvas widthsCanvas("c_residual_widths",
                       "Residual sigma and RMS by true position", 1000, 700);
  residualWidths.Draw("A");
  residualWidths.GetXaxis()->SetLimits(kSplineMinimumTrackY,
                                       kSplineMaximumTrackY);
  residualWidths.GetYaxis()->SetTitle("Width [#mum]");
  residualWidths.SetMinimum(0.0);
  residualWidths.SetMaximum(600.0);
  std::unique_ptr<TF1> linearResidualResolutionFit;
  std::array<std::unique_ptr<TF1>, 3> resolutionComponents;
  std::unique_ptr<TPaveText> resolutionFitParameters;
  TGraphErrors resolutionFitData;
  TGraph resolutionFitCurve;
  std::array<TGraph, 3> resolutionFitComponentCurves;
  TMultiGraph resolutionFitPlot;
  std::unique_ptr<TCanvas> resolutionFitCanvas;
  std::unique_ptr<TLegend> resolutionFitLegend;
  std::unique_ptr<TPaveText> resolutionFitUnavailableText;
  resolutionFitData.SetName("linear_residual_sigma_fit_data");
  resolutionFitData.SetTitle(
      "Linear-fit residual width;True track y [cm];Width [#mum]");
  resolutionFitData.SetMarkerStyle(22);
  resolutionFitData.SetMarkerColor(kMagenta + 1);
  resolutionFitData.SetLineColor(kMagenta + 1);
  for (int point = 0; point < linearResidualSigmaGraph.GetN(); ++point) {
    double x = 0.0;
    double y = 0.0;
    linearResidualSigmaGraph.GetPoint(point, x, y);
    const int fitDataPoint = resolutionFitData.GetN();
    resolutionFitData.SetPoint(fitDataPoint, x, y);
    resolutionFitData.SetPointError(
        fitDataPoint, linearResidualSigmaGraph.GetErrorX(point),
        linearResidualSigmaGraph.GetErrorY(point));
  }
  resolutionFitPlot.SetName("linear_residual_resolution_fit_plot");
  resolutionFitPlot.SetTitle(
      "Linear-fit residual resolution components;"
      "True track y [cm];Width [#mum]");
  resolutionFitPlot.Add(&resolutionFitData, "P");
  if (linearResidualSigmaGraph.GetN() >= 3) {
    double minimumWidth = std::numeric_limits<double>::infinity();
    double maximumWidth = 0.0;
    for (int point = 0; point < linearResidualSigmaGraph.GetN(); ++point) {
      double x = 0.0;
      double width = 0.0;
      linearResidualSigmaGraph.GetPoint(point, x, width);
      minimumWidth = std::min(minimumWidth, width);
      maximumWidth = std::max(maximumWidth, width);
    }

    linearResidualResolutionFit = std::make_unique<TF1>(
        "linear_residual_resolution_components_fit",
        residualResolutionComponents, kResolutionFitMinimumTrackY,
        kSplineMaximumTrackY, 5);
    linearResidualResolutionFit->SetTitle(
        "Linear-fit residual resolution components");
    linearResidualResolutionFit->SetParNames("n", "d", "c", "f", "k");
    constexpr double kPositiveParameterMinimum = 1.0e-6;
    const double maximumAbsTrackY =
        std::abs(kResolutionFitMinimumTrackY);
    linearResidualResolutionFit->SetParameters(
        40.0, maximumWidth * maximumWidth / maximumAbsTrackY, 70.0, 100.0,
        20.0);
    linearResidualResolutionFit->SetParLimits(0, 10.000001, 100.0);
    linearResidualResolutionFit->SetParLimits(
        1, 20,
        100.0 * maximumWidth * maximumWidth / maximumAbsTrackY);
    linearResidualResolutionFit->SetParLimits(
        2, kPositiveParameterMinimum,
        std::max(10.0 * maximumWidth, 140.0));
    linearResidualResolutionFit->SetParLimits(
        3, kPositiveParameterMinimum, std::max(10.0 * maximumWidth, 200.0));
    linearResidualResolutionFit->SetParLimits(4, 1.0, 100.0);
    const int fitStatus =
        linearResidualSigmaGraph.Fit(linearResidualResolutionFit.get(),
                                     "QRS0");
    if (fitStatus == 0) {
      linearResidualResolutionFit->SetLineColor(kOrange + 7);
      linearResidualResolutionFit->SetLineWidth(3);
      linearResidualResolutionFit->Draw("SAME");
      constexpr std::array<const char*, 3> componentNames = {
          "primary_ionization_resolution", "diffusion_resolution",
          "electronics_resolution"};
      constexpr std::array<int, 3> componentColors = {
          kBlue + 1, kGreen + 2, kCyan + 2};
      const std::array<double, 5> componentParameters = {
          linearResidualResolutionFit->GetParameter(0),
          linearResidualResolutionFit->GetParameter(1),
          linearResidualResolutionFit->GetParameter(2),
          linearResidualResolutionFit->GetParameter(3),
          linearResidualResolutionFit->GetParameter(4)};
      const std::array<double (*)(double*, double*), 3> componentFunctions = {
          primaryIonizationResolution, diffusionResolution,
          electronicsResolution};
      for (std::size_t component = 0;
           component < resolutionComponents.size(); ++component) {
        resolutionComponents[component] = std::make_unique<TF1>(
            componentNames[component], componentFunctions[component],
            kResolutionFitMinimumTrackY, kSplineMaximumTrackY,
            component == 2 ? 3 : 1);
        if (component == 2) {
          resolutionComponents[component]->SetParameters(
              componentParameters[2], componentParameters[3],
              componentParameters[4]);
        } else {
          resolutionComponents[component]->SetParameter(
              0, componentParameters[component]);
        }
        resolutionComponents[component]->SetLineColor(
            componentColors[component]);
        resolutionComponents[component]->SetLineStyle(2);
        resolutionComponents[component]->SetLineWidth(2);
        resolutionComponents[component]->Draw("SAME");
      }

      resolutionFitCurve.SetName("linear_residual_resolution_fit_curve");
      resolutionFitCurve.SetLineColor(kOrange + 7);
      resolutionFitCurve.SetLineWidth(3);
      constexpr std::array<const char*, 3> componentCurveNames = {
          "primary_ionization_resolution_curve",
          "diffusion_resolution_curve", "electronics_resolution_curve"};
      for (std::size_t component = 0;
           component < resolutionFitComponentCurves.size(); ++component) {
        resolutionFitComponentCurves[component].SetName(
            componentCurveNames[component]);
        resolutionFitComponentCurves[component].SetLineColor(
            componentColors[component]);
        resolutionFitComponentCurves[component].SetLineStyle(2);
        resolutionFitComponentCurves[component].SetLineWidth(2);
      }
      constexpr int kResolutionCurveSamples = 200;
      for (int sample = 0; sample < kResolutionCurveSamples; ++sample) {
        const double x =
            kResolutionFitMinimumTrackY +
            (kSplineMaximumTrackY - kResolutionFitMinimumTrackY) * sample /
                (kResolutionCurveSamples - 1);
        resolutionFitCurve.SetPoint(
            sample, x, linearResidualResolutionFit->Eval(x));
        for (std::size_t component = 0;
             component < resolutionComponents.size(); ++component) {
          resolutionFitComponentCurves[component].SetPoint(
              sample, x, resolutionComponents[component]->Eval(x));
        }
      }
      resolutionFitPlot.Add(&resolutionFitCurve, "L");
      for (auto& componentCurve : resolutionFitComponentCurves) {
        resolutionFitPlot.Add(&componentCurve, "L");
      }

      constexpr std::array<const char*, 5> parameterNames = {
          "n", "d", "c", "f", "k"};
      resolutionFitParameters =
          std::make_unique<TPaveText>(0.14, 0.62, 0.43, 0.89, "NDC");
      resolutionFitParameters->SetFillStyle(0);
      resolutionFitParameters->SetBorderSize(0);
      resolutionFitParameters->SetTextAlign(12);
      resolutionFitParameters->SetTextSize(0.03);
      const std::streamsize previousPrecision = std::cout.precision();
      std::cout << "Linear-fit residual resolution model parameters:"
                << std::setprecision(6);
      for (std::size_t parameter = 0; parameter < parameterNames.size();
           ++parameter) {
        const double value =
            linearResidualResolutionFit->GetParameter(parameter);
        const double error =
            linearResidualResolutionFit->GetParError(parameter);
        std::ostringstream label;
        label << parameterNames[parameter] << " = " << std::setprecision(4)
              << value << " #pm " << error;
        resolutionFitParameters->AddText(label.str().c_str());
        std::cout << " " << parameterNames[parameter] << " = " << value
                  << " +/- " << error;
      }
      std::cout.precision(previousPrecision);
      std::cout << " (fit uncertainties)." << std::endl;
      resolutionFitParameters->Draw();
    } else {
      std::cerr << "Resolution-components fit to linear-fit residual widths "
                   "failed (status "
                << fitStatus << ")." << std::endl;
      linearResidualResolutionFit.reset();
    }
  } else {
    std::cerr << "At least three linear-fit residual width points are "
                 "required for the resolution-components fit; found "
              << linearResidualSigmaGraph.GetN() << "." << std::endl;
  }
  resolutionFitCanvas = std::make_unique<TCanvas>(
      "c_resolution_components_fit",
      "Linear-fit residual resolution components", 1000, 700);
  resolutionFitCanvas->cd();
  resolutionFitPlot.Draw("A");
  if (resolutionFitData.GetN() > 0) {
    resolutionFitPlot.GetXaxis()->SetLimits(kResolutionFitMinimumTrackY,
                                            kSplineMaximumTrackY);
    resolutionFitPlot.GetYaxis()->SetTitle("Width [#mum]");
    resolutionFitPlot.SetMinimum(0.0);
    resolutionFitPlot.SetMaximum(600.0);
  }
  resolutionFitLegend = std::make_unique<TLegend>(0.62, 0.73, 0.89, 0.89);
  resolutionFitLegend->AddEntry(&resolutionFitData, "Linear residuals", "p");
  if (linearResidualResolutionFit) {
    resolutionFitLegend->AddEntry(&resolutionFitCurve, "Total fit", "l");
    resolutionFitLegend->AddEntry(
        &resolutionFitComponentCurves[0],
        "Primary ionization: n e^{-n|x|}", "l");
    resolutionFitLegend->AddEntry(
        &resolutionFitComponentCurves[1], "Diffusion: #sqrt{d|x|}", "l");
    resolutionFitLegend->AddEntry(
        &resolutionFitComponentCurves[2], "Electronics: c+fe^{kx}", "l");
  } else {
    resolutionFitUnavailableText =
        std::make_unique<TPaveText>(0.2, 0.5, 0.8, 0.6, "NDC");
    resolutionFitUnavailableText->SetFillStyle(0);
    resolutionFitUnavailableText->SetBorderSize(0);
    resolutionFitUnavailableText->SetTextAlign(22);
    resolutionFitUnavailableText->AddText(
        "Resolution fit unavailable; see macro output for details");
    resolutionFitUnavailableText->Draw();
  }
  resolutionFitLegend->Draw();
  resolutionFitCanvas->Modified();
  resolutionFitCanvas->Update();
  widthsCanvas.cd();

  TLegend widthsLegend(0.62, 0.76, 0.88, 0.89);
  widthsLegend.AddEntry(&residualSigmaGraph, "Spline Crystal Ball #sigma",
                        "p");
  widthsLegend.AddEntry(&linearResidualSigmaGraph,
                        "Linear Crystal Ball #sigma", "p");
  widthsLegend.AddEntry(&residualRmsGraph, "Residual RMS", "lp");
  if (linearResidualResolutionFit) {
    widthsLegend.AddEntry(linearResidualResolutionFit.get(),
                          "#sqrt{(n e^{-n|x|})^{2}+d|x|+(c+fe^{kx})^{2}} fit",
                          "l");
    widthsLegend.AddEntry(resolutionComponents[0].get(),
                          "Primary ionization: n e^{-n|x|}", "l");
    widthsLegend.AddEntry(resolutionComponents[1].get(),
                          "Diffusion: #sqrt{d|x|}", "l");
    widthsLegend.AddEntry(resolutionComponents[2].get(),
                          "Electronics: c+fe^{kx}", "l");
  }
  widthsLegend.Draw();
  widthsCanvas.Modified();
  widthsCanvas.Update();
  outputFile->cd();
  if (residualSigmaGraph.Write() <= 0 ||
      linearResidualSigmaGraph.Write() <= 0 ||
      residualRmsGraph.Write() <= 0 || residualWidths.Write() <= 0 ||
      (linearResidualResolutionFit &&
       linearResidualResolutionFit->Write() <= 0) ||
      (linearResidualResolutionFit &&
       (resolutionComponents[0]->Write() <= 0 ||
        resolutionComponents[1]->Write() <= 0 ||
        resolutionComponents[2]->Write() <= 0)) ||
      (resolutionFitParameters && resolutionFitParameters->Write() <= 0) ||
      (resolutionFitCanvas &&
       (resolutionFitData.Write() <= 0 || resolutionFitPlot.Write() <= 0 ||
        (linearResidualResolutionFit &&
         (resolutionFitCurve.Write() <= 0 ||
          resolutionFitComponentCurves[0].Write() <= 0 ||
          resolutionFitComponentCurves[1].Write() <= 0 ||
          resolutionFitComponentCurves[2].Write() <= 0)) ||
        resolutionFitCanvas->Write() <= 0)) ||
      widthsCanvas.Write() <= 0) {
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
  std::ifstream driftVelocityPlotFile(driftVelocityPlotPath,
                                      std::ios::binary);
  if (!driftVelocityPlotFile) {
    std::cerr << "Failed to create drift-velocity plot: "
              << driftVelocityPlotPath << std::endl;
    return;
  }
  std::cout << "Spatial-resolution plot saved to " << plotPath << std::endl;
  std::cout << "Crystal Ball time-resolution plot saved to "
            << timeResolutionPlotPath << " (" << timeResolutionGraph.GetN()
            << " successful slice fits)." << std::endl;
  std::cout << "MPV graph saved to " << rootPath << std::endl;
  std::cout << "Residual distributions saved to " << residualPlotPath
            << " (" << reconstructedEntries << " events; " << rejectedEntries
            << " events could not be inverted uniquely with the spline; "
            << linearReconstructedEntries << " events reconstructed with "
            << "the linear fit; " << linearRejectedEntries
            << " linear reconstructions fell outside the fit range)."
            << std::endl;
  std::cout << "Residual sigma and RMS plot saved to " << widthsPlotPath
            << std::endl;
  std::cout << "Resolution-fit canvas c_resolution_components_fit saved to "
            << rootPath << std::endl;
  std::cout << "Crystal Ball residual MPV canvas and graphs saved to "
            << rootPath << std::endl;
  std::cout << "Drift-velocity plot saved to " << driftVelocityPlotPath
            << " (spline derivative and linear-fit slope converted to "
               "cm/us; "
            << garfieldCurveCount << " Garfield gas-mixture curves)."
            << std::endl;
}
