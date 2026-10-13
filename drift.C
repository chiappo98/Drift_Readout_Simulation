#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TFile.h>
#include <ROOT/RDataFrame.hxx>

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

#include <nlohmann/json.hpp>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"
// #include "Garfield/ComponentAnsys123.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh" 
#include "Garfield/ViewField.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/ViewMedium.hh"

using namespace Garfield;

struct SimulationConfig {
  int shaperOrder;
  double shapingTimeNs;
  double gain;
  double threshold;
  double vSense;
  double vField;
  double vMylar;
  double cellW;
  double cellY;
  double dSense;
  double dField;
  double argonPercent;
  double co2Percent;
};

struct SignalInductionConfig {
  bool enabled;
  int numberOfEvents;
  double yMin;
  double yMax;
  bool printCumulative;
  bool plotSignals;
};

struct SimulationFlags {
  bool plotDrift;
  bool plotField;
  bool plotGasPar;
  bool simTrack;
  bool simDriftRKF;
  SignalInductionConfig signalInduction;
};

double shapingTimeNs = 25.0;

SimulationConfig readConfig(const std::string& configPath) {
  std::ifstream input(configPath);
  if (!input) {
    throw std::runtime_error("Cannot open configuration file: " + configPath);
  }

  const auto config = nlohmann::json::parse(input);
  if (!config.is_object() || !config.contains("shaper") ||
      !config["shaper"].is_number_integer() ||
      !config.contains("shaping_time_ns") ||
      !config["shaping_time_ns"].is_number() ||
      !config.contains("gain_mV/fC") || !config["gain_mV/fC"].is_number() ||
      !config.contains("threshold_fC") || !config["threshold_fC"].is_number() ||
      !config.contains("vSense") || !config["vSense"].is_number() ||
      !config.contains("vField") || !config["vField"].is_number() ||
      !config.contains("vMylar") || !config["vMylar"].is_number() ||
      !config.contains("cell_w") || !config["cell_w"].is_number() ||
      !config.contains("cell_y") || !config["cell_y"].is_number() ||
      !config.contains("dSense") || !config["dSense"].is_number() ||
      !config.contains("dField") || !config["dField"].is_number() ||
      !config.contains("argon_percent") || !config["argon_percent"].is_number() ||
      !config.contains("co2_percent") || !config["co2_percent"].is_number()) {
    throw std::runtime_error(
        "Configuration is missing a required simulation parameter or has an invalid type.");
  }

  const SimulationConfig result{
      config["shaper"].get<int>(),
      config["shaping_time_ns"].get<double>(),
      config["gain_mV/fC"].get<double>(),
      config["threshold_fC"].get<double>(),
      config["vSense"].get<double>(),
      config["vField"].get<double>(),
      config["vMylar"].get<double>(),
      config["cell_w"].get<double>(),
      config["cell_y"].get<double>(),
      config["dSense"].get<double>(),
      config["dField"].get<double>(),
      config["argon_percent"].get<double>(),
      config["co2_percent"].get<double>()};
  if ((result.shaperOrder != 1 && result.shaperOrder != 3) ||
      result.shapingTimeNs < 25.0 || result.shapingTimeNs > 250.0 ||
      result.gain <= 0.0 ||
      result.cellW <= 0.0 || result.cellY <= 0.0 ||
      result.dSense <= 0.0 || result.dField <= 0.0 ||
      result.argonPercent < 0.0 || result.argonPercent > 100.0 ||
      result.co2Percent < 0.0 || result.co2Percent > 100.0 ||
      std::abs(result.argonPercent + result.co2Percent - 100.0) > 1.e-6) {
    throw std::runtime_error(
        "Invalid configuration: check shaper, shaping time, event count, position, positive geometry values, and gas percentages summing to 100.");
  }
  return result;
}

SimulationFlags readSimulationFlags(const std::string& configPath) {
  std::ifstream input(configPath);
  if (!input) {
    throw std::runtime_error("Cannot open simulation configuration file: " + configPath);
  }

  const auto config = nlohmann::json::parse(input);
  if (!config.is_object() || !config.contains("plotDrift") ||
      !config["plotDrift"].is_boolean() || !config.contains("plotField") ||
      !config["plotField"].is_boolean() || !config.contains("plotGasPar") ||
      !config["plotGasPar"].is_boolean() || !config.contains("simTrack") ||
      !config["simTrack"].is_boolean() || !config.contains("simDriftRKF") ||
      !config["simDriftRKF"].is_boolean() || !config.contains("signalInduction") ||
      !config["signalInduction"].is_object() ||
      !config["signalInduction"].contains("enabled") ||
      !config["signalInduction"]["enabled"].is_boolean() ||
      !config["signalInduction"].contains("number_of_events") ||
      !config["signalInduction"]["number_of_events"].is_number_integer() ||
      !config["signalInduction"].contains("y_min") ||
      !config["signalInduction"]["y_min"].is_number() ||
      !config["signalInduction"].contains("y_max") ||
      !config["signalInduction"]["y_max"].is_number() ||
      !config["signalInduction"].contains("print_cumulative") ||
      !config["signalInduction"]["print_cumulative"].is_boolean() ||
      !config["signalInduction"].contains("plot_signals") ||
      !config["signalInduction"]["plot_signals"].is_boolean()) {
    throw std::runtime_error(
        "Simulation config must define the five mode flags and all signalInduction settings with the expected types.");
  }

  const SignalInductionConfig signalInduction{
      config["signalInduction"]["enabled"].get<bool>(),
      config["signalInduction"]["number_of_events"].get<int>(),
      config["signalInduction"]["y_min"].get<double>(),
      config["signalInduction"]["y_max"].get<double>(),
      config["signalInduction"]["print_cumulative"].get<bool>(),
      config["signalInduction"]["plot_signals"].get<bool>()};
    if (signalInduction.numberOfEvents < 1 || signalInduction.yMin < 0.0 ||
      signalInduction.yMin > 1.0 || signalInduction.yMax < 0.0 ||
      signalInduction.yMax > 1.0 || signalInduction.yMin > signalInduction.yMax) {
    throw std::runtime_error(
      "signalInduction.number_of_events must be >= 1 and y_min/y_max must satisfy 0 <= y_min <= y_max <= 1.");
  }

  return {config["plotDrift"].get<bool>(),
          config["plotField"].get<bool>(),
          config["plotGasPar"].get<bool>(),
          config["simTrack"].get<bool>(),
          config["simDriftRKF"].get<bool>(),
          signalInduction};
}

// bool readTransferFunction(Sensor& sensor) {
//   std::ifstream infile;
//   infile.open("mdt_elx_delta.txt", std::ios::in);
//   if (!infile) {
//     std::cerr << "Could not read delta response function.\n";
//     return false;
//   }
//   std::vector<double> times;
//   std::vector<double> values;
//   while (!infile.eof()) {
//     double t = 0., f = 0.;
//     infile >> t >> f;
//     if (infile.eof() || infile.fail()) break;
//     times.push_back(1.e3 * t);
//     values.push_back(f);
//   }
//   infile.close();
//   sensor.SetTransferFunction(times, values);
//   return true;
// }

// This function simulates the response of the front-end electronics (shaper)
double transfer1ord(double time) {
  return (time / shapingTimeNs) * exp(1.0 - time / shapingTimeNs);
}

double transfer3ord(double time) {
  return pow(time / shapingTimeNs, 3) * exp(3 * (1.0 - time / shapingTimeNs));
}

double compute_mean(std::vector<double> v) {
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  return sum / v.size();
}

double compute_std(std::vector<double> v) {
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  double mean = sum / v.size();
  double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
  return std::sqrt(sq_sum / v.size() - mean * mean);
}

int main(int argc, char* argv[]) {
  if (argc > 3) {
    std::cerr << "Usage: " << argv[0]
              << " [config/config.json] [config/simulation.json]" << std::endl;
    return 1;
  }
  const std::string configPath = argc >= 2 ? argv[1] : "config/config.json";
  const std::string simulationConfigPath = argc == 3 ? argv[2] : "config/simulation.json";
  SimulationConfig config;
  SimulationFlags simulationFlags;
  try {
    config = readConfig(configPath);
    simulationFlags = readSimulationFlags(simulationConfigPath);
  } catch (const std::exception& error) {
    std::cerr << "Configuration error: " << error.what() << std::endl;
    return 1;
  }
  shapingTimeNs = config.shapingTimeNs;

  int rootArgc = 1;
  char* rootArgv[] = {argv[0], nullptr};
  TApplication app("app", &rootArgc, rootArgv);

  // Setup the gas.
  MediumMagboltz gas;
  gas.SetTemperature(293.15); //298.15
  gas.SetPressure(763.7); // 5 mBar overpressure
  gas.SetComposition("ar", config.argonPercent, "co2", config.co2Percent);
  gas.Initialise(true);
  gas.EnablePenningTransfer();
  // Set the Penning transfer efficiency.
  // constexpr double rPenning = 0.51;
  // constexpr double lambdaPenning = 0.;
  // gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");

  MediumMagboltz gas1;
  gas1.SetTemperature(293.15); //298.15
  gas1.SetPressure(763.7); // 5 mBar overpressure
  gas1.SetComposition("ar", 70., "co2", 30.);
  gas1.Initialise(true);
  gas1.EnablePenningTransfer();

  MediumMagboltz gas2;
  gas2.SetTemperature(293.15); //298.15
  gas2.SetPressure(763.7); // 5 mBar overpressure
  gas2.SetComposition("ar", 85., "co2", 15.);
  gas2.Initialise(true);
  gas2.EnablePenningTransfer();

  // Load gas file with diffusion information (pre-calculated with Magboltz)
  if (!gas.LoadGasFile("/home/student/Software/GARFIELDPP_with_GEANT4_from_source_20260510/Examples/DriftDune/gas_files/ar_co2_90_10.gas")) {
      std::cerr << "Error: Could not load gas file. Ensure it exists." << std::endl;
      return 1;
  }
  if (!gas1.LoadGasFile("/home/student/Software/GARFIELDPP_with_GEANT4_from_source_20260510/Examples/DriftDune/gas_files/ar_co2_70_30.gas")) {
      std::cerr << "Error: Could not load gas file. Ensure it exists." << std::endl;
      return 1;
  }
  if (!gas2.LoadGasFile("/home/student/Software/GARFIELDPP_with_GEANT4_from_source_20260510/Examples/DriftDune/gas_files/ar_co2_85_15.gas")) {
    std::cerr << "Error: Could not load gas file. Ensure it exists." << std::endl;
    return 1;
  }
  // Load the ion mobilities.
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
  gas.LoadIonMobility("IonMobility_CO2+_CO2.txt");

  // Make a component with analytic electric field.
  ComponentAnalyticField cmp;
  cmp.SetMedium(&gas);
  // Wire radius [cm]
  const double dSense = config.dSense;
  const double dField = config.dField;
  // Cell width
  const double cell_w = config.cellW; // y direction
  // Cell thickness
  const double cell_t = config.cellY; // x direction
  // Voltages
  const double vSense = config.vSense;
  const double vField = config.vField;
  const double vMylar = config.vMylar;
  // Add the wire in the centre.
  cmp.AddWire(0, 0, dSense, vSense, "s");
  cmp.AddWire(0, -cell_w/2, dField, vField);
  cmp.AddWire(0, cell_w/2, dField, vField);
  // Add the Mylar planes.
  cmp.AddPlaneX(-cell_t/2, vMylar);
  cmp.AddPlaneX(cell_t/2, vMylar);

  // cmp.EnablePeriodicityY();
  cmp.SetPeriodicityY(cell_w);

  // Create the sensor.
  Sensor sensor(&cmp);
  sensor.AddElectrode(&cmp, "s");
  sensor.SetArea(-cell_t/2,-cell_w/2, -cell_w, cell_t/2, cell_w/2, cell_w);

  ViewDrift driftView;
  ViewField fieldView;
  ViewCell cellView;
  ViewMedium mediumView, mediumView1, mediumView2;


  // --------------------- Simulations ----------------------------------------------------------------------------------//
  const bool plotDrift = simulationFlags.plotDrift;
  const bool plotField = simulationFlags.plotField;
  const bool plotGasPar = simulationFlags.plotGasPar;
  const bool simTrack = simulationFlags.simTrack;
  const bool simDriftRKF = simulationFlags.simDriftRKF;
  const bool signalInduction = simulationFlags.signalInduction.enabled;
  const bool plotSignals = simulationFlags.signalInduction.plotSignals;
  // -------------------------------------------------------------------------------------------------------------------//

  fieldView.SetSensor(&sensor);
  cellView.SetComponent(&cmp);
  mediumView.SetMedium(&gas);
  mediumView1.SetMedium(&gas1);
  mediumView2.SetMedium(&gas2);
  
  if (plotField) {
    // TCanvas* cf = new TCanvas("cf", "multi-cell view", 300, 1000);
    // // fieldView.SetArea(-cell_t/2,-cell_w, 0, cell_t/2, cell_w, 0);
    // fieldView.SetCanvas(cf);
    // fieldView.Plot("e","COLZ");

    // Plot Electric Field Map
    TCanvas* cf1 = new TCanvas("cf1", "single cell E-field view", 700, 800);
    fieldView.EnableAutoRange();
    // fieldView.SetArea(-cell_t/2,-cell_w/2, cell_t/2, cell_w/2);
    fieldView.SetCanvas(cf1);
    fieldView.Plot("e","COLZ");
    gPad->SetLogz(1);

    // Plot Electric Potential Map
    TCanvas* cf11 = new TCanvas("cf11", "single cell E-field view", 700, 800);
    fieldView.EnableAutoRange();
    // fieldView.SetArea(-cell_t/2,-cell_w/2, cell_t/2, cell_w/2);
    fieldView.SetCanvas(cf11);
    fieldView.Plot("p","COLZ");
    gPad->SetLogz(0);


    TCanvas* cf2 = new TCanvas("cf2", "E-field lines", 700, 800);
    fieldView.SetCanvas(cf2);
    cellView.SetCanvas(cf2);
    std::vector<double> line_x, line_y, line_z;
    fieldView.EqualFluxIntervals(
      -0.99*cell_t/2, -0.99*cell_w/2, 0.0,
      -0.99*cell_t/2, +0.99*cell_w/2, 0.0,
      line_x, line_y, line_z, 50
    );
    fieldView.EqualFluxIntervals(
      +0.99*cell_t/2, -0.99*cell_w/2, 0.0,
      +0.99*cell_t/2, +0.99*cell_w/2, 0.0,
      line_x, line_y, line_z, 50
    );
    fieldView.EqualFluxIntervals(
      -0.99*cell_t/2, -0.99*cell_w/2, 0.0,
      +0.99*cell_t/2, -0.99*cell_w/2, 0.0,
      line_x, line_y, line_z, (int)(50 * cell_t / cell_w)
    );
    fieldView.EqualFluxIntervals(
      -0.99*cell_t/2, +0.99*cell_w/2, 0.0,
      +0.99*cell_t/2, +0.99*cell_w/2, 0.0,
      line_x, line_y, line_z, (int)(50 * cell_t / cell_w)
    );
    fieldView.PlotFieldLines(line_x, line_y, line_z);
    cellView.Plot2d();


    // Generate Drift Velocity Map
    Medium *med;
    const double map_spacing = 0.005;
    const int np_x = (int)(cell_t / map_spacing); 
    const int np_y = (int)(cell_w / map_spacing);
    double xp = 0., yp = 0., zp = 0., t_fin = 0.;
    double e_gain = 0., e_loss = 0., drift_t = 0., i_drift_t = 0., dt_sigma = 0.,
           i_dt_sigma = 0., path_len = 0., av_ne = 0., av_ni = 0.;
    int drift_st = 0, i_drift_st = 0;
    std::array<double, 3> e_field = {0., 0., 0.};
    std::array<double, 3> b_field = {0., 0., 0.};
    std::array<double, 3> drift_v = {0., 0., 0.};
    double v_mag = 0.;

    TH2D* h_v = new TH2D("h_v", "Drift Velocity", np_x, -cell_t/2, cell_t/2, np_y, -cell_w/2, cell_w/2);

      // map the electron response
    for (int cx = 0; cx <= np_x; cx++) {
      xp = -cell_t/2 + cx * map_spacing;
      for (int cy = 0; cy <= np_y; cy++) {
        yp = -cell_w/2 + cy * map_spacing;

        // 7-1 get the E and B field components for
        sensor.MagneticField(xp, yp, zp, b_field[0], b_field[1], b_field[2], drift_st);
        sensor.ElectricField(xp, yp, zp, e_field[0], e_field[1], e_field[2], med, drift_st);
        if (drift_st != 0) {
          // drift_map.FillLUT(blank_vec); // fill with the blank vector
          continue;
        }
        med->ElectronVelocity(e_field[0], e_field[1], e_field[2], b_field[0],
                              b_field[1], b_field[2], drift_v[0], drift_v[1],
                              drift_v[2]);
        v_mag = 1e4*std::sqrt(drift_v[0]*drift_v[0]+drift_v[1]*drift_v[1]+drift_v[2]*drift_v[2]);
        h_v->Fill(xp, yp, v_mag);
      }
    }

    gPad->SetLogz(0);
    TCanvas* c_v = new TCanvas("c_v", "Drift Velocity Map", 700, 800);
    c_v->cd();
    c_v->SetRightMargin(0.15);
    h_v->SetStats(false);
    h_v->Draw("colz");
    h_v->SetMaximum(50.);
    h_v->GetXaxis()->SetTitle("x (cm)");
    h_v->GetYaxis()->SetTitle("z (cm)");
    c_v->SetRealAspectRatio();

    // Output File Path
    const std::string outFilePath = "../plot/mylar/0.7w-0.6h_1550m-1950f_noComsol.root";
    // const std::string outFilePath = "../plot/mesh/0.7w-0.6h_1550m-1950f_0.2p-0.05d.root";

    // Save to file
    TFile f( outFilePath.c_str(), "RECREATE");
    h_v->Write();
    std::cout << "Saved plot to File " << outFilePath << std::endl; 


  }
  if (plotGasPar) {
    TCanvas* cm = new TCanvas("cm", "Electron Drift Velocity", 1800, 600);
    cm->Divide(3,1);

    TPad* p1 = (TPad*)cm->cd(1);
    mediumView.SetCanvas(p1);
    mediumView.PlotElectronVelocity('e');
    TLatex l1; l1.SetTextAlign(13); l1.DrawLatexNDC(0.3, 0.5, "Ar-CO_{2} 90-10"); 
    TPad* p2 = (TPad*)cm->cd(2);
    mediumView2.SetCanvas(p2);
    mediumView2.PlotElectronVelocity('e');
    TLatex l2; l2.SetTextAlign(13); l2.DrawLatexNDC(0.3, 0.5, "Ar-CO_{2} 85-15"); 
    TPad* p3 = (TPad*)cm->cd(3);
    mediumView1.SetCanvas(p3);
    mediumView1.PlotElectronVelocity('e');
    TLatex l3; l3.SetTextAlign(13); l3.DrawLatexNDC(0.3, 0.5, "Ar-CO_{2} 70-30"); 

    TCanvas* cA = new TCanvas("cA", "Transport Parameters", 1800, 600);
    cA->Divide(3,1);
    
    TPad* pA1 = (TPad*)cA->cd(1);
    mediumView.SetCanvas(pA1);
    mediumView.PlotElectronDiffusion('e');
    TLatex lA1; lA1.SetTextAlign(13); lA1.DrawLatexNDC(0.3, 0.3, "Ar-CO_{2} 90-10"); 
    TPad* pA2 = (TPad*)cA->cd(2);
    mediumView2.SetCanvas(pA2);
    mediumView2.PlotElectronDiffusion('e');
    TLatex lA2; lA2.SetTextAlign(13); lA2.DrawLatexNDC(0.3, 0.3, "Ar-CO_{2} 85-15"); 
    TPad* pA3 = (TPad*)cA->cd(3);
    mediumView1.SetCanvas(pA3);
    mediumView1.PlotElectronDiffusion('e');
    TLatex lA3; lA3.SetTextAlign(13); lA3.DrawLatexNDC(0.3, 0.3, "Ar-CO_{2} 70-30"); 
  }
  if (plotDrift) {

    TCanvas* c1 = new TCanvas("c1", "Wire Amplification", 800, 600);

    ViewDrift* driftView_wire = new ViewDrift();
    driftView_wire->SetCanvas(c1);

    AvalancheMicroscopic* aval_wire = new AvalancheMicroscopic();
    aval_wire->SetSensor(&sensor);
    aval_wire->EnablePlotting(driftView_wire, 1);

    // Start electron close to the wire (0.02 cm) to see the avalanche quickly
    double x0 = -0.7*cell_t/2, y0 = -cell_w/2, z0 = 0.0, t0 = 0.0, e0 = 0.0;
    aval_wire->AvalancheElectron(x0, y0, z0, t0, e0, 0, 0, 0);

    int ne, ni;
    aval_wire->GetAvalancheSize(ne, ni);
    std::cout << "Number of final electrons created near the wire: " << ne << std::endl;

    if(ne > 0) driftView_wire->Plot(true);
    c1->SaveAs("AvalancheMicroscopic_Wire_Line.png");

  }
  if (simTrack) {

    TrackHeed *track = new TrackHeed();
    track->SetParticle("electron");
    track->SetEnergy(450e6); // eV
    track->SetSensor(&sensor);

    // Initial position and direction
    double x0 = -0.7*cell_t/2, y0 = -cell_w/2, z0 = 0.0, t0 = 0.0;
    double dx0 = 0.0, dy0 = 1.0, dz0 = 0.0;

    // TH1I *h_Ncl = new TH1I("h_Ncl_A", "Number of primary clusters", 100, 0, 100);
    TH1D *h_clX = new TH1D("h_clX_A", "Cluster position X; X [cm]", 100, -2.0, 2.0);
    TH1D *h_clY = new TH1D("h_clY_A", "Cluster position Y; Y [cm]", 100, -2.0, 2.0);
    TH1D *h_clZ = new TH1D("h_clZ_A", "Cluster position Z; Z [cm]", 100, -2.0, 2.0);
    // TH1D *h_clDistance = new TH1D("h_clDistance_A", "Cluster distance; Distance [cm]", 100, 0.0, 0.5);

    TH1I *h_Ncl = new TH1I("h_Ncl", "Number of primary clusters; Number of primary cluster", 100, 0, 100);
    TH1D *h_clDistance = new TH1D("h_clDistance", "Cluster distance; Distance [cm]", 100, 0.0, 0.5);

    TH1I *h_Nsec = new TH1I("h_Nsec", "Number of secondary electrons; Number of secondary electrons", 40, 0, 40);
    TH1D *h_dE = new TH1D("h_dE", "Deposited energy in primary cluster; Energy deposit in primary clusters [eV]", 250, 0.0, 1000.0);
    TH1D *h_tot_dE = new TH1D("h_tot_dE", "Total deposited energy by particle; Total Energy deposit per particle [eV]", 200, 0.0, 10000.0);
    TH2D *h_dE_vs_Nsec = new TH2D("h_dE_vs_Nsec", "dE vs N secondary; Energy [eV]; N Secondary", 100, 0.0, 1000.0, 40, 0, 40);


    int Ntimes = 10000;
    for(int iSim = 0; iSim < Ntimes; iSim++) {
        if(iSim % 1000 == 0) std::cout << " Simulating event " << iSim << " / " << Ntimes << std::endl;
        
        track->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
        const auto& clusters = track->GetClusters();
        
        h_Ncl->Fill(clusters.size());
        double totE = 0.0;

        for(size_t iCl = 1; iCl < clusters.size(); iCl++) {
            h_clX->Fill(clusters[iCl].x);
            h_clY->Fill(clusters[iCl].y);
            h_clZ->Fill(clusters[iCl].z);

            double dX = clusters[iCl].x - clusters[iCl-1].x;
            double dY = clusters[iCl].y - clusters[iCl-1].y;
            double dZ = clusters[iCl].z - clusters[iCl-1].z;
            h_clDistance->Fill(sqrt(dX*dX + dY*dY + dZ*dZ));

            int nElectrons = clusters[iCl].electrons.size();
            double energy = clusters[iCl].energy;
            h_Nsec->Fill(nElectrons);
            h_dE->Fill(energy);
            h_dE_vs_Nsec->Fill(energy, nElectrons);
            totE += energy;
        }
        h_tot_dE->Fill(totE);
    }

    TCanvas *cA = new TCanvas("cA", "Particle Info", 1200, 500);

    cA->Divide(3, 1);
    cA->cd(1); h_clX->Draw(); 
    cA->cd(2); h_clY->Draw();
    cA->cd(3); h_clZ->Draw(); 
    cA->SaveAs("particle_info.png");


    TCanvas *cC = new TCanvas("cC", "Part C: Energy", 1600, 800);
    cC->Divide(3, 2);
    cC->cd(1); h_Ncl->Draw(); TF1 *f_Ncl = new TF1("f_Ncl", "gaus"); h_Ncl->Fit(f_Ncl, "Q"); h_Ncl->Draw();
    cC->cd(2); gPad->SetLogy(); h_Nsec->Draw();
    cC->cd(3); h_clDistance->Draw(); TF1 *f_clDistance = new TF1("f_clDistance", "expo"); h_clDistance->Fit(f_clDistance, "Q");
    cC->cd(4); h_dE->GetXaxis()->SetRangeUser(0.0, 100.0); h_dE->Draw();
    double E_mpv = h_tot_dE->GetBinCenter(h_tot_dE->GetMaximumBin());
    TF1 *f_E = new TF1("f_E", "gaus", 0.6*E_mpv, 1.4*E_mpv);
    cC->cd(5); h_tot_dE->Fit(f_E, "RQ"); h_tot_dE->Draw();
    cC->cd(6); gPad->SetLogz(); h_dE_vs_Nsec->Draw("colz");
    
    cC->SaveAs("track_info.png");

    double cl_size_mu, cl_size_std, E_mu, E_std;
    cl_size_mu = f_Ncl->GetParameter(1);
    cl_size_std = f_Ncl->GetParameter(2);
    E_mu = f_E->GetParameter(1);
    E_std = f_E->GetParameter(2);

    std::cout << "ratio $\sigma / \mu$ for: " << std::endl;
    std::cout << "cluster size -> " << cl_size_std/cl_size_mu << std::endl;
    std::cout << "total energy -> " << E_std/E_mu << std::endl;
  }
  if (simDriftRKF) {
    int N_events = 1;
    int repetitions = 2;

    // Define ionizing track (muon)
    TrackHeed* track_rkf = new TrackHeed();
    track_rkf->SetParticle("muon");
    track_rkf->SetEnergy(170e9);
    track_rkf->SetSensor(&sensor);

    // Define Drift method (RKF is deterministic/smooth)
    DriftLineRKF* drift_RKF = new DriftLineRKF();
    drift_RKF->SetSensor(&sensor);
    drift_RKF->EnableIonTail();

    ViewDrift* driftView_rkf = new ViewDrift();
    
    if (N_events == 1) {
         
      // Connect plotting to both track and drift objects
      track_rkf->EnablePlotting(driftView_rkf);
      drift_RKF->EnablePlotting(driftView_rkf);

      repetitions = 1;
    }
         
    TGraphErrors *gr = new TGraphErrors();
    
    for (double e = 0.; e < N_events; ++e ) {
      // Initial track parameters: passing horizontally above the wire
      if (N_events==1) e=.9;
      double x1 = -cell_t/2, y1 = -(e/N_events) * cell_w/2, z1 = 0.0, t1 = 0.0;
      double dx1 = 1.0, dy1 = 0.0, dz1 = 0.0;
      // double x1 = -(e/N_events) * cell_t/2, y1 = -cell_w/2, z1 = 0.0, t1 = 0.0;
      // double dx1 = 0.0, dy1 = 1.0, dz1 = 0.0;
  
      std::vector<double> t_spread, distance;
      for (int r = 0; r < repetitions; ++r) {
        track_rkf->NewTrack(x1, y1, z1, t1, dx1, dy1, dz1);

        // Extract and drift electrons from primary clusters
        const auto& clusters_rkf = track_rkf->GetClusters();
        for (const auto& cl : clusters_rkf) {
            for (const auto& electron : cl.electrons) {
                drift_RKF->DriftElectron(electron.x, electron.y, electron.z, electron.t);
                // drift_RKF->DriftIon(electron.x, electron.y, electron.z, electron.t);
            }
        }
        t_spread.push_back(drift_RKF->GetArrivalTimeSpread()); // ns
        distance.push_back((e/N_events) * cell_w/2);
      }
      if (N_events>1) {
        gr->SetPoint((int)(e), compute_mean(distance), compute_mean(t_spread));
        gr->SetPointError((int)(e), compute_std(distance), compute_std(t_spread));
      }
      if (N_events == 1) {
        TCanvas* cRKF = new TCanvas("cRKF", "RKF Drift Lines", 600, 600);
        driftView_rkf->SetCanvas(cRKF);

        const bool twod = true;
        driftView_rkf->SetArea(-cell_t/2,-cell_w/2, cell_t/2, cell_w/2);
        driftView_rkf->Plot(twod);
        cellView.SetCanvas(cRKF);
        cellView.Plot2d();
        cRKF->SaveAs("RKF_Drift_Line.png");

        std::cout << " " << std::endl;
        std::cout << "Arrival time spread (ns): " << drift_RKF->GetArrivalTimeSpread() << std::endl;
        std::cout << "Gain: " << drift_RKF->GetGain() << std::endl;
        std::cout << "Loss: " << drift_RKF->GetLoss() << std::endl;

      }
    }
    if (N_events>1) {
      TCanvas *cts = new TCanvas("cts", "Diffusion", 700, 500);
      cts->cd();
      gr->GetXaxis()->SetTitle("distance (cm)");
      gr->GetYaxis()->SetTitle("time spread (ns)");
      gr->Draw("APL");
    }
  }
  if (signalInduction) {

    std::vector<double> v_totalCharge, v_q_peak, v_t_peak, v_t_threshold;
    std::vector<double> v_x1, v_y1, v_z1, v_dx1, v_dy1, v_dz1;

    double nbins = 8000;
    double tStep = 0.5;
    double tStart = 0.0;

    TCanvas* c1 = nullptr;
    TH1* ht = nullptr;

    sensor.SetTimeWindow(tStart, tStep, nbins);

    // Define ionizing track (muon)
    TrackHeed* track_ind = new TrackHeed();
    track_ind->SetParticle("muon");
    track_ind->SetEnergy(170e9);
    track_ind->SetSensor(&sensor);

    // Define Drift method (RKF is deterministic/smooth)
    DriftLineRKF* drift_ind = new DriftLineRKF();
    drift_ind->SetSensor(&sensor);
    drift_ind->EnableIonTail();
    drift_ind->UseWeightingPotential();
    
    ViewDrift* driftView_rkf = new ViewDrift();
    
    // Connect plotting to both track and drift objects
    track_ind->EnablePlotting(driftView_rkf);
    drift_ind->EnablePlotting(driftView_rkf);

    std::random_device randomSeed;
    std::mt19937 randomEngine(randomSeed());
    std::uniform_real_distribution<double> positionDistribution(
      simulationFlags.signalInduction.yMin, simulationFlags.signalInduction.yMax);

    for (int iEvent = 0; iEvent < simulationFlags.signalInduction.numberOfEvents; ++iEvent) {
      const int progressStep = std::max(1, simulationFlags.signalInduction.numberOfEvents / 10);
      if (iEvent % progressStep == 0) {
        std::cout << "Simulating event " << iEvent + 1 << " / "
                  << simulationFlags.signalInduction.numberOfEvents << std::endl;
      }
      if (iEvent > 0) sensor.ClearSignal();
      const double position = simulationFlags.signalInduction.yMin == simulationFlags.signalInduction.yMax
                  ? simulationFlags.signalInduction.yMin
                  : positionDistribution(randomEngine);
      const double x1 = -cell_t/2;
      const double y1 = - position * cell_w/2;
      const double z1 = 0.0;
      const double t1 = 0.0;
      const double dx1 = 1.0;
      const double dy1 = 0.0;
      const double dz1 = 0.0;
      track_ind->NewTrack(x1, y1, z1, t1, dx1, dy1, dz1);
      const auto& clusters_ind = track_ind->GetClusters();
      for (const auto& cl : clusters_ind) {
        for (const auto& electron : cl.electrons) {
          drift_ind->DriftElectron(electron.x, electron.y, electron.z, electron.t);
          drift_ind->DriftIon(electron.x, electron.y, electron.z, electron.t);
        }
      }


      if (plotSignals) {
        c1 = new TCanvas("c1", "Induction Signals", 2000, 500);
        c1->Divide(4, 1);

        c1->cd(1);
        driftView_rkf->SetCanvas((TPad*)c1->GetPad(1));
        const bool twod = true;
        driftView_rkf->SetArea(-cell_t/2,-cell_w/2, cell_t/2, cell_w/2);
        driftView_rkf->Plot(twod);
        cellView.SetCanvas((TPad*)c1->GetPad(1));
        cellView.Plot2d();

        TLatex* tex = new TLatex();
        tex->SetNDC(true);
        tex->DrawLatex(0.2, 0.85, Form("Arrival time spread: %.1f ns", drift_ind->GetArrivalTimeSpread()));

        c1->cd(2);
        gStyle->SetLineWidth(1);
        sensor.PlotSignal("s",(TPad*)c1->GetPad(2),"tei");

        TLegend* leg = new TLegend(0.65, 0.15, 0.8, 0.33);
        TPad* pad2 = (TPad*)c1->GetPad(2);
        ht = (TH1*)pad2->GetPrimitive("hSignal__0");
        TH1* he = (TH1*)pad2->GetPrimitive("hSignalElectrons__0");
        TH1* hi = (TH1*)pad2->GetPrimitive("hSignalIons__0");
        ht->SetLineWidth(2);
        he->SetLineWidth(2);
        hi->SetLineWidth(2);
        leg->AddEntry(ht, "total", "l");
        leg->AddEntry(he, "electron", "l");
        leg->AddEntry(hi, "ion", "l");
        leg->Draw();
        pad2->Modified();
        pad2->Update();
      }

      // --- COMPUTE CUMULATIVE CHARGE DISTRIBUTION
      TH1* hch = (TH1*) ht->Clone("hch");
      std::vector<double> totCharge = {0};
      double totalCharge = 0;
      for (int i = 1; i <= hch->GetNbinsX(); ++i) {
        totCharge.push_back(totCharge[i-1] - ht->GetBinContent(i)*tStep);
        totalCharge -= ht->GetBinContent(i)*tStep;
        hch -> SetBinContent(i, totCharge[i]);
      }
      if (simulationFlags.signalInduction.printCumulative) {
        double charge99 = 0.9 * totalCharge;
        auto it = std::find_if(totCharge.begin(), totCharge.end(),
                            [charge99](double v) { return v > charge99; });
        int bin99 = std::distance(totCharge.begin(), it);

        TCanvas* c2 = new TCanvas("c2", "Cumulative Signal Charge", 600, 500);
        c2->cd();
        hch->Draw();

        double ychmin = gPad->GetUymin(); // current y-range of the pad (after Draw())
        double ychmax = gPad->GetUymax();
        double t99 = hch->GetBinCenter(bin99);
        TLine* line = new TLine(t99, ychmin, t99, ychmax);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->SetLineStyle(2); // dashed
        line->Draw();
        gPad->Modified();
        gPad->Update();

        std::cout << " " << std::endl;
        std::cout << "t99: " << t99 << std::endl;
        std::cout << "charge99: " << charge99 << std::endl;
      }

      sensor.SetTransferFunction(config.shaperOrder == 1 ? transfer1ord : transfer3ord);
      sensor.ConvoluteSignals();
      double eventPeak = 0.0;
      int eventPeakBin = 0;
      int eventThresholdBin = -1;
      for (int i = 0; i < nbins; ++i) {
        const double signal = sensor.GetSignal("s", i);
        if (i == 0 || signal < eventPeak) {
          eventPeak = signal;
          eventPeakBin = i;
        }
        if (eventThresholdBin < 0 && signal < config.threshold) {
          eventThresholdBin = i;
        }
      }

      v_totalCharge.push_back(totalCharge);
      v_q_peak.push_back(eventPeak);
      v_t_peak.push_back(tStart + (eventPeakBin + 0.5) * tStep);
      v_t_threshold.push_back(eventThresholdBin < 0
                                    ? -1.0
                                    : tStart + (eventThresholdBin + 0.5) * tStep);
      v_x1.push_back(x1);
      v_y1.push_back(y1);
      v_z1.push_back(z1);
      v_dx1.push_back(dx1);
      v_dy1.push_back(dy1);
      v_dz1.push_back(dz1);



      // --- 3. PROCESSED SIGNAL (Electronics Transfer Function) ---
      const double gain = config.gain;
      const double threshold = config.threshold;
    
      TH1* ht3 = nullptr;
      TPad* pad3 = nullptr;
      if (plotSignals) {
        c1->cd(3);
        sensor.PlotSignal("s",(TPad*)c1->GetPad(3));
        pad3 = (TPad*)c1->GetPad(3);
        ht3 = (TH1*)pad3->GetPrimitive("hSignal__0");
        ht3->SetLineWidth(2);

        double ymin = ht3->GetMinimum();
        double ymax = ht3->GetMaximum();
        pad3->Update();
        TGaxis* axis_mV = new TGaxis(
            pad3->GetUxmax(), pad3->GetUymin(),
            pad3->GetUxmax(), pad3->GetUymax(),
            ymin * gain, ymax * gain,
            510, "+L");
        axis_mV->SetTitle("Signal (mV)");
        axis_mV->SetTitleOffset(1.3);
        axis_mV->SetLineColor(kRed);
        axis_mV->SetLabelColor(kRed);
        axis_mV->SetTitleColor(kRed);
        axis_mV->Draw();

      }

      // FIND peak TIME AND CHARGE and time-at-threshold
  
      int bin_below = -1;
      for (int i = 1; i <= ht3->GetNbinsX(); ++i) {
          if (ht3->GetBinContent(i) < threshold) {
              bin_below = i;
              break;
          }
      }

      if (plotSignals) {
        TLatex* tex1 = new TLatex();
        tex1->SetNDC(true);
        std::vector<std::string> lines = {
            Form("Q signal: %.1f fC", totalCharge),
            Form("Q measured (peak): %.1f fC", eventPeak),
            Form("t@peak: %.1f ns", v_t_peak[-1]),
            Form("t@threshold: %.1f ns", v_t_threshold[-1])
        };
        double y0 = 0.3;
        double dy = 0.05;
        for (size_t i = 0; i < lines.size(); ++i) {
          tex1->DrawLatex(0.4, y0 - i*dy, lines[i].c_str());
        }

        pad3->SetTicks(1, 0);
        pad3->SetRightMargin(0.15);
        pad3->SetBorderMode(0);
        pad3->Modified();
        pad3->Update();
      }

      if (plotSignals) {
        if (ht->Integral() != 0) ht->Scale(1.0 / ht->Integral());
        if (ht3->Integral() != 0) ht3->Scale(3.0 / ht3->Integral());

        c1->cd(4);
        TPad* pad4 = (TPad*)c1->GetPad(4);
        ht->SetLineColor(kBlue);
        ht->SetLineWidth(2);
        ht3->SetLineColor(kRed);
        ht3->SetLineWidth(2);

        double y4max = std::max(ht->GetMaximum(), ht3->GetMaximum());
        ht->SetMaximum(y4max * 1.1);
        ht->Draw("HIST");
        ht3->Draw("HIST SAME");

        TLegend* leg4 = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg4->AddEntry(ht, "current signal", "l");
        leg4->AddEntry(ht3, "shaper signal", "l");
        leg4->Draw();

        pad4->Modified();
        pad4->Update();
      }

    }

    // dopo il loop, crea l'RDataFrame dai vettori
    ROOT::RDataFrame df(v_totalCharge.size());
    auto df2 = df.Define("totalCharge", [&](ULong64_t i){ return v_totalCharge[i]; }, {"rdfentry_"})
                .Define("q_peak",      [&](ULong64_t i){ return v_q_peak[i]; }, {"rdfentry_"})
                .Define("t_peak",      [&](ULong64_t i){ return v_t_peak[i]; }, {"rdfentry_"})
                .Define("t_threshold", [&](ULong64_t i){ return v_t_threshold[i]; }, {"rdfentry_"})
                .Define("x1", [&](ULong64_t i){ return v_x1[i]; }, {"rdfentry_"})
                .Define("y1", [&](ULong64_t i){ return v_y1[i]; }, {"rdfentry_"})
                .Define("z1", [&](ULong64_t i){ return v_z1[i]; }, {"rdfentry_"})
                .Define("dx1", [&](ULong64_t i){ return v_dx1[i]; }, {"rdfentry_"})
                .Define("dy1", [&](ULong64_t i){ return v_dy1[i]; }, {"rdfentry_"})
                .Define("dz1", [&](ULong64_t i){ return v_dz1[i]; }, {"rdfentry_"});
    df2.Snapshot("results", "../results/results.root"); // salva su file come TTree, riapribile con RDataFrame


    sensor.ClearSignal();
  }


  app.Run();
}
