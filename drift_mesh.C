#include <TApplication.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TFile.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"

// #include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ComponentComsol.hh"

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

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  // Setup the gas.
  MediumMagboltz gas;
  gas.SetTemperature(293.15); //298.15
  gas.SetPressure(763.7); // 5 mBar overpressure
  gas.SetComposition("ar", 90., "co2", 10.);
  gas.Initialise(true);
  gas.EnablePenningTransfer();
  // Set the Penning transfer efficiency.
  // constexpr double rPenning = 0.51;
  // constexpr double lambdaPenning = 0.;
  // gas.EnablePenningTransfer(rPenning, lambdaPenning, "ar");

  // Load gas file with diffusion information (pre-calculated with Magboltz)
  if (!gas.LoadGasFile("/home/student/Software/GARFIELDPP_with_GEANT4_from_source_20260510/Examples/DriftDune/gas_files/ar_co2_90_10.gas")) {
      std::cerr << "Error: Could not load gas file. Ensure it exists." << std::endl;
      return 1;
  }

  // Load the ion mobilities.
  gas.LoadIonMobility("IonMobility_Ar+_Ar.txt");
  gas.LoadIonMobility("IonMobility_CO2+_CO2.txt");

  // Output File Path
  // const std::string outFilePath = "../plot/mylar/0.7w-0.6h_1550m-1950f.root";
  const std::string outFilePath = "../plot/mesh/0.7w-0.6h_1550m-1950f_0.2p-0.05d.root";


  // Path to get the FEM map
  // const std::string mapMeshPath = "../maps/mylar/0.7w-0.6h_1550m-1950f/";
  const std::string mapMeshPath = "../maps/mesh/0.7w-0.6h_1550m-1950f_0.2p-0.05d/";

  // Load the field map for woven mesh
  ComponentComsol fmMesh;
  fmMesh.Initialise(mapMeshPath + "mesh.mphtxt", mapMeshPath + "dielectrics.dat",
                    mapMeshPath + "field.txt", "m");
  
  fmMesh.EnableMirrorPeriodicityX();
  fmMesh.EnableMirrorPeriodicityY();
  fmMesh.PrintRange();

  int nMaterials = fmMesh.GetNumberOfMaterials();
  std::cout << "GetNumberOfMaterials = "<< nMaterials << std::endl;
  for (unsigned int i = 0; i < nMaterials; ++i) {
    const double eps = fmMesh.GetPermittivity(i);
    std::cout << "eps = "<< eps << std::endl;
    if(eps==1) fmMesh.SetMedium(i, &gas);
  }
  fmMesh.PrintMaterials();

  ComponentComsol fm;
  fm.SetImportMaterial(5e10); // Only importe domain that has the gas medium
  fm.Initialise(mapMeshPath + "mesh.mphtxt", mapMeshPath + "dielectrics.dat",
                mapMeshPath + "field.txt", "m");

  // Wire radius [cm]
  const double dSense = 0.002;
  const double dField = 0.01;
  // Cell width
  const double cell_w = 1.4; // y direction
  // Cell thickness
  const double cell_t = 1.2; // x direction
  // Voltages
  const double vSense = 0.;
  const double vField = -1950.;
  const double vMylar = -1500.;


  e

  // Create the sensor.
  Sensor sensor(&fmMesh);
  sensor.AddElectrode(&fm, "s");
  // sensor.SetArea(-cell_t/2,-cell_w/2, -cell_w, cell_t/2, cell_w/2, cell_w);

  ViewDrift driftView;
  ViewField fieldView;
  ViewCell cellView;
  ViewMedium mediumView;

  constexpr bool plotField = true;

  fieldView.SetSensor(&sensor);
  // cellView.SetComponent(&cmp);
  mediumView.SetMedium(&gas);
  

  if (plotField) {

    fieldView.SetPlane(0, -1, 0, 0, 0.35, 0);
    TCanvas* cf1 = new TCanvas("cf1", "single cell E-field view", 700, 600);
    // fieldView.EnableAutoRange();
    fieldView.SetArea(-cell_w/2,-cell_t/2, cell_w/2, cell_t/2);
    fieldView.SetCanvas(cf1);
    fieldView.Plot("e","COLZ");
    gPad->SetLogz(1);

    fieldView.SetPlane(0, -1, 0, 0, 0.3, 0);
    TCanvas* cf = new TCanvas("cf", "single cell Electric Potential view", 700, 600);
    fieldView.SetCanvas(cf);
    fieldView.Plot("p","zcol");

    // TCanvas* cf2 = new TCanvas("cf2", "E-field lines", 700, 800);
    // fieldView.SetCanvas(cf2);
    // cellView.SetCanvas(cf2);
    // std::vector<double> line_x, line_y, line_z;
    // fieldView.EqualFluxIntervals(
    //   -0.99*cell_t/2, -0.99*cell_w/2, 0.0,
    //   -0.99*cell_t/2, +0.99*cell_w/2, 0.0,
    //   line_x, line_y, line_z, 50
    // );
    // fieldView.EqualFluxIntervals(
    //   +0.99*cell_t/2, -0.99*cell_w/2, 0.0,
    //   +0.99*cell_t/2, +0.99*cell_w/2, 0.0,
    //   line_x, line_y, line_z, 50
    // );
    // fieldView.EqualFluxIntervals(
    //   -0.99*cell_t/2, -0.99*cell_w/2, 0.0,
    //   +0.99*cell_t/2, -0.99*cell_w/2, 0.0,
    //   line_x, line_y, line_z, (int)(50 * cell_t / cell_w)
    // );
    // fieldView.EqualFluxIntervals(
    //   -0.99*cell_t/2, +0.99*cell_w/2, 0.0,
    //   +0.99*cell_t/2, +0.99*cell_w/2, 0.0,
    //   line_x, line_y, line_z, (int)(50 * cell_t / cell_w)
    // );
    // fieldView.PlotFieldLines(line_x, line_y, line_z);
    // cellView.Plot2d();

    // //Generate Drift Velocity Map
    const double map_spacing = 0.005;
    const int np_z = (int)(cell_t / map_spacing); 
    const int np_x = (int)(cell_w / map_spacing);
    double xp = 0., yp = 0., zp = 0., t_fin = 0.;
    double e_gain = 0., e_loss = 0., drift_t = 0., i_drift_t = 0., dt_sigma = 0.,
           i_dt_sigma = 0., path_len = 0., av_ne = 0., av_ni = 0.;
    int drift_st = 0, i_drift_st = 0;
    std::array<double, 3> e_field = {0., 0., 0.};
    std::array<double, 3> b_field = {0., 0., 0.};
    std::array<double, 3> drift_v = {0., 0., 0.};
    double v_mag = 0.;

    TH2D* h_v = new TH2D("h_v", "Drift Velocity", np_x, -cell_w/2, cell_w/2, np_z, -cell_t/2, cell_t/2);

      // map the electron response
    for (int cx = 0; cx <= np_x; cx++) {
      xp = -cell_w/2 + cx * map_spacing;
      for (int cz = 0; cz <= np_z; cz++) {
        zp = -cell_t/2 + cz * map_spacing;

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
        h_v->Fill(xp, zp, v_mag);
      }
    }

    gPad->SetLogz(0);
    TCanvas* c_v = new TCanvas("c_v", "Drift Velocity Map", 700, 600);
    c_v->cd();
    c_v->SetRightMargin(0.15);
    h_v->SetStats(false);
    h_v->Draw("colz");
    h_v->SetMaximum(50.);
    h_v->GetXaxis()->SetTitle("x (cm)");
    h_v->GetYaxis()->SetTitle("z (cm)");
    c_v->SetRealAspectRatio();


    // Save to file
    TFile f( outFilePath.c_str(), "RECREATE");
    h_v->Write();
    std::cout << "Saved plot to File " << outFilePath << std::endl; 
  
  }



  app.Run();
}
