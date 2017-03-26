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

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"

int main( int argc, char** argv ) {

    // retrieve simulations and data
    std::string path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumData.root";
    TFile fileData(path.c_str(), "READ");
    if (!fileData.IsOpen()) { std::cout << "Zombie fileData!\n"; return -1; }

    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_2nbb.root";
    TFile fileSim2nbb(path.c_str(), "READ"); 
    if (!fileSim2nbb.IsOpen()) { std::cout << "Zombie fileSim2nbb!\n"; return -1; }

    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_2nbbLV.root";
    TFile fileSim2nbbLV(path.c_str(), "READ");
    if (!fileSim2nbbLV.IsOpen()) { std::cout << "Zombie fileSim2nbbLV!\n"; return -1; }

    // remember: pointers filled with TDirectoryFile::Get are non-owning!
    TH1D* hDataBEGe;   fileData.GetObject("energyBEGeAll", hDataBEGe);
    TH1D* hDataCOAX;   fileData.GetObject("energyEnrCoaxAll", hDataCOAX);
    TH1F* h2nbbBEGe;   fileSim2nbb.GetObject("energy_BEGe", h2nbbBEGe);
    TH1F* h2nbbCOAX;   fileSim2nbb.GetObject("energy_COAX", h2nbbCOAX);
    TH1F* h2nbbLVBEGe; fileSim2nbbLV.GetObject("energy_BEGe", h2nbbLVBEGe);
    TH1F* h2nbbLVCOAX; fileSim2nbbLV.GetObject("energy_COAX", h2nbbLVCOAX);

    if (!hDataBEGe or !hDataCOAX or !h2nbbBEGe or
        !h2nbbCOAX or !h2nbbLVBEGe or !h2nbbLVCOAX ) { std::cout << "There's at least one zombie TH1D!\n"; return -1; }
    
    // create binning
    std::vector<int> ubin(7500);
    int k = 1;
    std::generate(ubin.begin(), ubin.end(), [&](){ return k++; } );

    // convert in std::vector
    // spectrum components
    const int nComp = 2;
    std::vector<int> vDataBEGe, vDataCOAX;
    std::vector<std::vector<double>> vSimCOAX(nComp), vSimBEGe(nComp);

    for ( int i = 0; i < 7500; ++i ) {
        vDataBEGe.push_back(hDataBEGe->GetBinContent(i+1));
        vDataCOAX.push_back(hDataCOAX->GetBinContent(i+1));
        
        vSimBEGe[0].push_back(h2nbbBEGe->GetBinContent(i+1));
        vSimCOAX[0].push_back(h2nbbCOAX->GetBinContent(i+1));
        vSimBEGe[1].push_back(h2nbbLVBEGe->GetBinContent(i+1));
        vSimCOAX[1].push_back(h2nbbLVCOAX->GetBinContent(i+1));
    }

    fileData.Close();
    fileSim2nbb.Close();
    fileSim2nbbLV.Close();

// ------------------------------------------------------------------------------
	// create Fit2nbbLV object
    Fit2nbbLV m("Fit2nbbLV");

    // set data
    m.SetBinning(ubin);
    m.SetDataBEGe(vDataBEGe);
    m.SetDataCOAX(vDataCOAX);
    m.SetSimBEGe(vSimBEGe);
    m.SetSimCOAX(vSimCOAX);

    m.SetFitRange(600, 1400);

    // set nicer style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
	BCLog::OpenLog("out/logBAT.txt", BCLog::detail, BCLog::detail);

	// set precision (number of samples in Markov chain)
	m.MCMCSetPrecision(BCEngineMCMC::kLow);
    
    // set parameter binning
    //m.SetNbins(1000);

	m.WriteMarkovChain(true);
	// normalize the posterior, i.e. integrate posterior over the full
	// parameter space
	m.SetIntegrationMethod(BCIntegrate::kIntGrid);
	m.Normalize();

    auto start = std::chrono::system_clock::now();
	// run MCMC and marginalize posterior w/r/t all parameters and all
	// combinations of two parameters
	m.MarginalizeAll();

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
	
    return 0;
}
