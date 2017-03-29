/* runfit.cxx
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 24/03/2017
 *
 */

#include "Fit2nbbLV.h"

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <memory>

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFrame.h"
#include "TLegend.h"

int main( int argc, char** argv ) {
    
/////////////////////////////////////////////
    int nBins = 750;
    int rangeUp = 2700;  // [keV]
    int rangeDown = 550; // [keV]
    BCEngineMCMC::Precision level(BCEngineMCMC::kLow);
/////////////////////////////////////////////

    TH1::AddDirectory(false);
    // retrieve data
    std::string path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumData.root";
    TFile fileData(path.c_str(), "READ");
    if (!fileData.IsOpen()) { std::cout << "Zombie fileData!\n"; return -1; }
    
    // vector holding all sim files
    std::vector<std::unique_ptr<TFile>> simFile;
    // [0] 2nbb
    // [1] 2nbbLV
    // [2] homLAr
    // [3] K40onFiberShroud
    
    // 2nbb
    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_2nbb.root";
    simFile.emplace_back( new TFile(path.c_str(), "READ") );

    // 2nbbLV
    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_2nbbLV.root";
    simFile.emplace_back( new TFile(path.c_str(), "READ") );

    // homLAr
    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_homLAr.root";
    simFile.emplace_back( new TFile(path.c_str(), "READ") );
    
    // K40onFiberShroud
    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_K40onFiberShroud.root";
    simFile.emplace_back( new TFile(path.c_str(), "READ") );
    
    // Bi212onFiberShroud
    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_Bi212onFiberShroud.root";
    simFile.emplace_back( new TFile(path.c_str(), "READ") );
    
    // Tl208onFiberShroud
    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_Tl208onFiberShroud.root";
    simFile.emplace_back( new TFile(path.c_str(), "READ") );

    // Pb214onFiberShroud
    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_Pb214onFiberShroud.root";
    simFile.emplace_back( new TFile(path.c_str(), "READ") );
    
    // Bi214onFiberShroud
    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_Bi214onFiberShroud.root";
    simFile.emplace_back( new TFile(path.c_str(), "READ") );
 
    for ( auto& f : simFile ) if (!f->IsOpen()) { std::cout << "At least one zombie simFile!\n"; return -1; }

    TH1D* hDataBEGe;   fileData.GetObject("energyBEGeAll", hDataBEGe);
    TH1D* hDataCOAX;   fileData.GetObject("energyEnrCoaxAll", hDataCOAX);
    
    std::vector<TH1F*> hSimBEGe;
    std::vector<TH1F*> hSimCOAX;
    
    for ( auto& f : simFile ) {
        hSimBEGe.push_back(dynamic_cast<TH1F*>(f->Get("energy_BEGe")));
        hSimCOAX.push_back(dynamic_cast<TH1F*>(f->Get("energy_COAX")));
    }

    fileData.Close();
    for ( auto& f : simFile ) f->Close();

    if (!hDataBEGe or !hDataCOAX) { std::cout << "There's at least one zombie data hist!\n"; return -1; }
    for ( auto& h : hSimBEGe ) if (!h) { std::cout << "There's at least one zombie simBEGe hist!\n"; return -1; }
    for ( auto& h : hSimCOAX ) if (!h) { std::cout << "There's at least one zombie simCOAX hist!\n"; return -1; }
    
    // create binning !!!NON-VARIABLE!!!
    std::vector<int> ubin(nBins);
    int k = 0;
    std::generate(ubin.begin(), ubin.end(), [&](){ return k+=7500/nBins; } );

    // convert in std::vector
    std::vector<int> vDataBEGe, vDataCOAX;
    std::vector<std::vector<double>> vSimCOAX(hSimCOAX.size()), vSimBEGe(hSimBEGe.size());

    // rebin if needed
    hDataBEGe->Rebin(7500/nBins);
    hDataCOAX->Rebin(7500/nBins);
    for ( unsigned int i = 0; i < hSimBEGe.size(); ++i ) {
        hSimBEGe[i]->Rebin(7500/nBins);
        hSimCOAX[i]->Rebin(7500/nBins);
    }

    for ( int i = 0; i < nBins; ++i ) {
        vDataBEGe.push_back(hDataBEGe->GetBinContent(i+1));
        vDataCOAX.push_back(hDataCOAX->GetBinContent(i+1));
        
        for ( unsigned int j = 0; j < hSimBEGe.size(); ++j ) {
            vSimBEGe[j].push_back(hSimBEGe[j]->GetBinContent(i+1));
            vSimCOAX[j].push_back(hSimCOAX[j]->GetBinContent(i+1));
        }
    }

// ------------------------------------------------------------------------------
	// create Fit2nbbLV object
    Fit2nbbLV m("Fit2nbbLV");

    // set data
    m.SetBinning(ubin);
    m.SetDataBEGe(vDataBEGe);
    m.SetDataCOAX(vDataCOAX);
    m.SetSimBEGe(vSimBEGe);
    m.SetSimCOAX(vSimCOAX);

    m.SetFitRange(rangeDown, rangeUp);

    // set nicer style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog("out/logBAT.txt", BCLog::detail, BCLog::detail);

	// set precision (number of samples in Markov chain)
	m.MCMCSetPrecision(level);
    
    // set parameter binning
    //m.SetNbins(1000);

	// normalize the posterior, i.e. integrate posterior over the full
	// parameter space
	m.SetIntegrationMethod(BCIntegrate::kIntDefault);
	m.Normalize();

    auto start = std::chrono::system_clock::now();
	// run MCMC and marginalize posterior w/r/t all parameters and all
	// combinations of two parameters
	m.MarginalizeAll();
	m.WriteMarkovChain(true);

	// run mode finding; by default using Minuit
	m.FindMode(m.GetBestFitParameters());
    
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start);

	// draw all marginalized distributions into a PDF file
	m.PrintAllMarginalized("out/Fit2nbbLV_plots.pdf");

	// create a new summary tool object, to print change from prior -> posterior
	BCSummaryTool summary(&m);
	summary.PrintKnowledgeUpdatePlots("out/Fit2nbbLV_update.pdf");

    // print all summary plots
	summary.PrintParameterPlot("out/Fit2nbbLV_parameters.pdf");
	summary.PrintCorrelationPlot("out/Fit2nbbLV_correlation.pdf");
	summary.PrintCorrelationMatrix("out/Fit2nbbLV_correlationMatrix.pdf");

	// calculate p-value
	// m -> CalculatePValue( m->GetBestFitParameters() );

	// print results of the analysis into a text file
	m.PrintResults("out/Fit2nbbLV_results.txt");

	BCLog::OutSummary("Exiting");

	// close log file
	BCLog::CloseLog();

    std::cout << "\n\ntime: " << elapsed.count() << "s\n";
// -----------------------------------------------------------------------------------

    // save results to draw them later
    std::vector<double> results = m.GetBestFitParameters();
    
    TFile outDrawFile("out/outHists.root", "RECREATE");
    
    hDataBEGe->SetName("hDataBEGe");
    hDataBEGe->Write();
    hDataCOAX->SetName("hDataCOAX");
    hDataCOAX->Write();

    // 2nbb
    hSimBEGe[0]->Scale(results[0]);
    hSimBEGe[0]->SetBins(nBins,0,7500);
    hSimBEGe[0]->SetName("h2nbbBEGe");

    hSimCOAX[0]->Scale(results[0]);
    hSimCOAX[0]->SetBins(nBins,0,7500);
    hSimCOAX[0]->SetName("h2nbbCOAX");
    
    // 2nbbLV
    hSimBEGe[1]->Scale(results[0]*results[1]*m.Getn2n1());
    hSimBEGe[1]->SetBins(nBins,0,7500);
    hSimBEGe[1]->SetName("h2nbbLVBEGe");
    
    hSimCOAX[1]->Scale(results[0]*results[1]*m.Getn2n1());
    hSimCOAX[1]->SetBins(nBins,0,7500);
    hSimCOAX[1]->SetName("h2nbbLVCOAX");
    
    // K42
    hSimBEGe[2]->Scale(results[2]);
    hSimBEGe[2]->SetName("hK42homLArBEGe");
    
    hSimCOAX[2]->Scale(results[2]);
    hSimCOAX[2]->SetName("hK42homLArCOAX");
    
    // K40
    hSimBEGe[3]->Scale(results[3]);
    hSimBEGe[3]->SetName("hK40onFiberShroudBEGe");
    
    hSimCOAX[3]->Scale(results[3]);
    hSimCOAX[3]->SetName("hK40onFiberShroudCOAX");
     
    // Bi212
    hSimBEGe[4]->Scale(results[4]);
    hSimBEGe[4]->SetName("hBi212onFiberShroudBEGe");
    
    hSimCOAX[4]->Scale(results[4]);
    hSimCOAX[4]->SetName("hBi212onFiberShroudCOAX");
    
    // Tl208
    hSimBEGe[5]->Scale(results[4]*m.GetBrRatioTl());
    hSimBEGe[5]->SetName("hTl208onFiberShroudBEGe");
    
    hSimCOAX[5]->Scale(results[4]*m.GetBrRatioTl());
    hSimCOAX[5]->SetName("hTl208onFiberShroudCOAX");
    
    // Pb214
    hSimBEGe[6]->Scale(results[5]);
    hSimBEGe[6]->SetName("hPb214onFiberShroudBEGe");
    
    hSimCOAX[6]->Scale(results[5]);
    hSimCOAX[6]->SetName("hPb214onFiberShroudCOAX");
    
    // Bi214
    hSimBEGe[7]->Scale(results[5]);
    hSimBEGe[7]->SetName("hBi214onFiberShroudBEGe");
    
    hSimCOAX[7]->Scale(results[5]);
    hSimCOAX[7]->SetName("hBi214onFiberShroudCOAX");
 
    for ( unsigned int i = 0; i < hSimBEGe.size(); ++i ) {
        hSimBEGe[i]->Write();
        hSimCOAX[i]->Write();
    }

// ----------------------------------------------------------------------------------
    
    auto draw = [&]( std::vector<TH1F*> v , TH1D* vd ,std::string type ) {
        
        TCanvas tmp(type.c_str(), type.c_str(), 1200, 700);
        TPad pad("pad", "pad", 0.0, 0.0, 1, 1);
        pad.SetMargin(0.07,0.05,0.1,0.05);
        tmp.cd();
        pad.Draw();
        pad.cd();
    
        vd->SetStats(false);
        vd->SetMarkerStyle(23);
        vd->GetXaxis()->SetRangeUser(rangeDown,rangeUp);
        vd->GetXaxis()->SetTitle("energy [keV]");
        vd->GetYaxis()->SetTitle("counts");
        vd->GetYaxis()->SetNdivisions(10);
        vd->Draw("P");

        v[0]->SetLineColor(kBlue);
        v[1]->SetLineColor(kBlue+2);  
        v[2]->SetLineColor(kGreen);
        v[3]->SetLineColor(kGreen+1);
        v[4]->SetLineColor(kGreen+2);
        v[5]->SetLineColor(kGreen+3);
        v[6]->SetLineColor(kGreen+4);
        v[7]->SetLineColor(kBlack);

        for ( auto& h : v ) h->Draw("SAME");
    
        std::string name = "sum" + type;
        TH1D sum(name.c_str(), name.c_str(), nBins, 0, 7500);
        for ( auto& h : v ) sum.Add(h);
        sum.SetLineColor(kRed);
        sum.Draw("SAME");

        TLegend leg(0.7,0.7,0.95,0.95);
        for ( auto& h : v ) leg.AddEntry(h,h->GetName(),"l");
        leg.AddEntry(&sum,sum.GetName(),"l");
        leg.Draw();
        
        pad.SetLogy();
        pad.SetGrid();
        name = "out/" + type + ".pdf";
        tmp.SaveAs(name.c_str());
    };

    draw( hSimBEGe, hDataBEGe, "BEGe" );
    draw( hSimCOAX, hDataCOAX, "COAX" );

    delete hDataBEGe;   
    delete hDataCOAX;   
    hSimBEGe.clear();
    hSimCOAX.clear();

    return 0;
}
