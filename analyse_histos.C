#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLatex.h>

#include <fstream>
#include <iostream>
#include <vector>


void analyse_histos() {

    // Cell width
    const double cell_w = 1.4; // y direction
    // Cell thickness
    const double cell_t = 1.2; // x direction
  
    // Input File Path
    const std::string mylarFilePath = "plot/mylar/0.7w-0.6h_1550m-1950f.root";
    const std::string meshFilePath = "plot/mesh/0.7w-0.6h_1550m-1950f_0.2p-0.05d.root";

    // Output File Path
    const std::string outFilePath = "plot/histo_analysis.root";

    TFile f1(mylarFilePath.c_str());
    TFile f2(meshFilePath.c_str());

    TH2D *h_v_mylar = (TH2D*)f1.Get("h_v");
    TH2D *h_v_mesh = (TH2D*)f2.Get("h_v");

    int nbins1x = h_v_mylar->GetNbinsX();
    int nbins1y = h_v_mylar->GetNbinsY();
    int nbins2x = h_v_mesh->GetNbinsX();
    int nbins2y = h_v_mesh->GetNbinsY();

    // Study drift velocity saturation

    std::vector<double> v_drift;
    int num = 0; 
    for (int x = 1; x <= nbins1x; ++x) {
        for (int y = 1; y <= nbins1y; ++y) {
            if ( h_v_mylar->GetBinContent(x, y) > 45. ){
                num += 1;
            }
        }
    }
    double sat_ratio_mylar = (double)num / (double)(nbins1x*nbins1y);

    num = 0; 
    for (int x = 1; x <= nbins2x; ++x) {
        for (int y = 1; y <= nbins2y; ++y) {
            if ( h_v_mesh->GetBinContent(x, y) > 45. ){
                num += 1;
            }
        }
    }
    double sat_ratio_mesh = (double)num / (double)(nbins2x*nbins2y);

    std::cout << "Saturated region (v > 45um/ns) - Mylar: " << sat_ratio_mylar << std::endl;
    std::cout << "Saturated region (v > 45um/ns) - Mesh: " << sat_ratio_mesh << std::endl;

    // Compare histograms

    if ( (nbins1x != nbins2x) | (nbins1y != nbins2y)) {
        std::cerr << "Mismatch in the bin number: " << nbins1x << " " << nbins1y << " - " << nbins2x << " " << nbins2y << std::endl;
        std::exit(1);
    }

    TH2D* h_v_comp = new TH2D("h_v_comp", "v mesh / v mylar", nbins1x, -cell_w/2, cell_w/2, nbins1y, -cell_t/2, cell_t/2);
    h_v_comp->SetDirectory(nullptr);
    
    for (int x = 1; x <= nbins1x; ++x) {
        for (int y = 1; y <= nbins1y; ++y) {
            double v1 = h_v_mylar->GetBinContent(x, y);
            double v2 = h_v_mesh->GetBinContent(x, y);
            if ( !std::isnan(v1) & !std::isnan(v2)){
                double diff = v2 / v1;
                if (std::fabs(diff) < 1e-5) diff = 0;
                h_v_comp->SetBinContent(x, y, diff);
            }
        }
    }

    TCanvas* c_v_comp = new TCanvas("c_v_comp", "Drift Velocity Comparison", 700, 600);
    c_v_comp->cd();
    h_v_comp->SetStats(false);
    h_v_comp->Draw("colz");
    h_v_comp->GetXaxis()->SetTitle("x (cm)");
    h_v_comp->GetYaxis()->SetTitle("z (cm)");



    // Save to file
    TFile f( outFilePath.c_str(), "RECREATE");
    h_v_comp->Write();
    std::cout << "Saved plot to File " << outFilePath << std::endl; 

}
