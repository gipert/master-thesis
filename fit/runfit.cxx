/* runfit.cxx
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 24/03/2017
 *
 * Variable binning:
 *   4 keV general binning, for the following gamma lines:
 *
 *   569    583    609    911    969    1001    1063
 *   1120   1173   1238   1332   1461   1525    1764
 *   2204   2614
 *
 *   we merge bins around +- 4keV ==> 1847 bins
 *
 */

#include "Fit2nbbLV.h"

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <memory>

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCModelOutput.h>

#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"

double GetPValue(Fit2nbbLV& m, bool save = false);

int main( int argc, char** argv ) {
    
/////////////////////////////////////////////
    const int rangeUp = 5300;  // [keV]
    const int rangeDown = 570; // [keV] above 39Ar Q-value
    BCEngineMCMC::Precision level(BCEngineMCMC::kMedium);
/////////////////////////////////////////////
    
    auto c_str = [](std::string s) { return s.c_str(); };

    TH1::AddDirectory(false);
    // retrieve data
    std::string path = std::string(std::getenv("GERDACPTDIR")) + "/data/";
    TFile fileData(c_str(path + "sumData.root"), "READ");
    if (!fileData.IsOpen()) { std::cout << "Zombie fileData!\n"; return -1; }
    
    // vector holding all sim files
    std::vector<std::unique_ptr<TFile>> simFile;
    // [0] 2nbb
    // [1] 2nbbLV
    // [2] homLAr
    // [3] K40onFiberShroud
    // [4] Bi212
    // [5] Tl208
    // [6] Pb214
    // [7] Bi214
    // [8] alphas
    // [9] nPlus
    // [10] pPlus
    // [11] Ac228holder
    // [12] Co60holder
    // [13] K40holder
    
    std::vector<std::string> simpath;

    simpath.push_back(path + "sumMaGe_2nbb.root");
    simpath.push_back(path + "sumMaGe_2nbbLV.root");
    simpath.push_back(path + "sumMaGe_K42homLAr.root");    
    simpath.push_back(path + "sumMaGe_K40fibers.root");    
    simpath.push_back(path + "sumMaGe_Bi212fibers.root");    
    simpath.push_back(path + "sumMaGe_Tl208fibers.root");
    simpath.push_back(path + "sumMaGe_Pb214fibers.root");   
    simpath.push_back(path + "sumMaGe_Bi214fibers.root");
    simpath.push_back(path + "alpha_model_run53-74.root");
    simpath.push_back(path + "sumMaGe_K42nPlus.root"); 
    simpath.push_back(path + "sumMaGe_K42pPlus.root");    
    simpath.push_back(path + "sumMaGe_Ac228holder.root");
    simpath.push_back(path + "sumMaGe_Co60holder.root"); 
    simpath.push_back(path + "sumMaGe_K40holder.root"); 
    simpath.push_back(path + "sumMaGe_Bi212holder.root"); 
    simpath.push_back(path + "sumMaGe_Tl208holder.root"); 
    simpath.push_back(path + "sumMaGe_Pb214holder.root"); 
    simpath.push_back(path + "sumMaGe_Bi214holder.root"); 
    simpath.push_back(path + "sumMaGe_K40cables.root"); 
    simpath.push_back(path + "sumMaGe_Bi212cables.root"); 
    simpath.push_back(path + "sumMaGe_Tl208cables.root"); 
    simpath.push_back(path + "sumMaGe_Pb214cables.root"); 
    simpath.push_back(path + "sumMaGe_Bi214cables.root"); 
    simpath.push_back(path + "sumMaGe_K40minishroud.root"); 
    simpath.push_back(path + "sumMaGe_Pa234minishroud.root"); 
    simpath.push_back(path + "sumMaGe_Bi207minishroud.root"); 
    simpath.push_back(path + "sumMaGe_Bi207cables.root"); 
    simpath.push_back(path + "sumMaGe_Bi207holder.root"); 

    // TODO: push_back new files here

    for ( const auto& p : simpath ) simFile.emplace_back( new TFile(p.c_str(), "READ") );
 
    for ( auto& f : simFile ) if (!f->IsOpen()) { std::cout << "At least one zombie simFile!\n"; return -1; }

    TH1D* hDataBEGetmp;   fileData.GetObject("energyBEGeAll", hDataBEGetmp);
    TH1D* hDataCOAXtmp;   fileData.GetObject("energyEnrCoaxAll", hDataCOAXtmp);
    
    std::vector<TH1*> hSimBEGetmp;
    std::vector<TH1*> hSimCOAXtmp;

    TH1* tmp;
    
    for ( auto& f : simFile ) {
        f->GetObject("energy_BEGe", tmp);
        hSimBEGetmp.push_back(tmp);
        f->GetObject("energy_COAX", tmp);
        hSimCOAXtmp.push_back(tmp);
    }
    
    fileData.Close();
    for ( auto& f : simFile ) f->Close();

    for ( auto& h : hSimBEGetmp ) if (!h) { std::cout << "There's at least one zombie simBEGe hist!\n"; return -1; }
    for ( auto& h : hSimCOAXtmp ) if (!h) { std::cout << "There's at least one zombie simCOAX hist!\n"; return -1; }
    
    // create !!!VARIABLE!!! binning
    
    const int nBins = 1847;
    std::vector<int> avoid = { 568, 572, 580, 584, 608, 612, 908, 
                               912, 968, 972, 1000, 1004, 1060, 1064, 
                               1120, 1172, 1176, 1236, 1240, 1332, 1460, 
                               1464, 1524, 1528, 1764, 2204, 2612, 2616 };

    std::vector<double> dbin(nBins+1);
    int k = 0, i = 0;
    while (1) { 
        if ( std::find( avoid.begin(), avoid.end(), k) == avoid.end() ) { 
            dbin[i] = k; 
            i++; 
        }
        k += 4;
        if ( k > 7500 ) break;
    }
    
    // my simulations are in MeV...
    hSimBEGetmp[0]->SetBins(7500,0,7500);
    hSimCOAXtmp[0]->SetBins(7500,0,7500);
    hSimBEGetmp[1]->SetBins(7500,0,7500);
    hSimCOAXtmp[1]->SetBins(7500,0,7500);
    
    // rebin the histograms
    // new histograms
    TH1D* hDataBEGe = dynamic_cast<TH1D*>(hDataBEGetmp->Rebin(nBins, "hDataBEGetmp", &dbin[0]));
    TH1D* hDataCOAX = dynamic_cast<TH1D*>(hDataCOAXtmp->Rebin(nBins, "hDataCOAXtmp", &dbin[0]));
    
    std::vector<TH1*> hSimBEGe;
    std::vector<TH1*> hSimCOAX;

    for ( unsigned int i = 0; i < hSimBEGetmp.size(); ++i ) {
        if ( i == 8 ) {
            hSimBEGe.push_back(dynamic_cast<TH1D*>(hSimBEGetmp[i]->Rebin(nBins, Form("hSimBEGe%i", i), &dbin[0])));
            hSimCOAX.push_back(dynamic_cast<TH1D*>(hSimCOAXtmp[i]->Rebin(nBins, Form("hSimCOAX%i", i), &dbin[0])));
        }
        else {
            hSimBEGe.push_back(dynamic_cast<TH1F*>(hSimBEGetmp[i]->Rebin(nBins, Form("hSimBEGe%i", i), &dbin[0])));
            hSimCOAX.push_back(dynamic_cast<TH1F*>(hSimCOAXtmp[i]->Rebin(nBins, Form("hSimCOAX%i", i), &dbin[0])));
        }
    }

    delete hDataBEGetmp;   
    delete hDataCOAXtmp;   
    hSimBEGetmp.clear();
    hSimCOAXtmp.clear();

    // convert in std::vector
    std::vector<int> vDataBEGe(nBins), vDataCOAX(nBins);
    std::vector<std::vector<double>> vSimCOAX(hSimCOAX.size()), vSimBEGe(hSimBEGe.size());

    // normalize alphas
    hSimBEGe[8]->Scale(1./hSimBEGe[8]->Integral());
    hSimCOAX[8]->Scale(1./hSimCOAX[8]->Integral());

    for ( int i = 0; i < nBins; ++i ) {
        vDataBEGe[i] = hDataBEGe->GetBinContent(i+1);
        vDataCOAX[i] = hDataCOAX->GetBinContent(i+1);
        
        for ( unsigned int j = 0; j < hSimBEGe.size(); ++j ) {
            vSimBEGe[j].push_back(hSimBEGe[j]->GetBinContent(i+1));
            vSimCOAX[j].push_back(hSimCOAX[j]->GetBinContent(i+1));
        }
    }

    delete hDataBEGe;   
    delete hDataCOAX;   
    hSimBEGe.clear();
    hSimCOAX.clear();

// ==========================================================================================================
	
    // create Fit2nbbLV object
    Fit2nbbLV model("Fit2nbbLV");
	// create a new summary tool object
	BCSummaryTool summary(&model);
    // create output class
    path = std::string(std::getenv("GERDACPTDIR")) + "/";
    BCModelOutput output(&model, c_str(path + "out/markowChains.root"));
    model.WriteMarkovChain(true);
    
    // set nicer style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
    BCLog::SetLogLevelFile(BCLog::detail);
    BCLog::SetLogLevelScreen(BCLog::summary);
	BCLog::OpenLog(c_str(path + "/out/logBAT.txt"));

	// set precision (number of samples in Markov chain)
	model.MCMCSetPrecision(level);
    
    // set data
    model.SetBinning(dbin);
    model.SetDataBEGe(vDataBEGe);
    model.SetDataCOAX(vDataCOAX);
    model.SetSimBEGe(vSimBEGe);
    model.SetSimCOAX(vSimCOAX);
    model.SetFitRange(rangeDown, rangeUp);

    // set parameter binning
    //m.SetNbins(1000);
    
	// run MCMC and marginalize posterior w/r/t all parameters and all
	// combinations of two parameters
    auto start = std::chrono::system_clock::now();
	model.MarginalizeAll();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start);	
    std::cout << "Summary : Time spent: " << elapsed.count() << "s.\n";

	// run mode finding, by default using Minuit
    BCLog::SetLogLevelScreen(BCLog::detail);
	model.FindMode(model.GetBestFitParameters());
    BCLog::SetLogLevelScreen(BCLog::summary);

    std::cout << std::endl;
    double pvalue = GetPValue(model);
    std::cout << "Summary : pValue = " << pvalue << std::endl;

    // OUTPUT
	// print results of the analysis into a text file
	model.PrintResults(c_str(path + "out/Fit2nbbLV_results.txt"));
	// draw all marginalized distributions into a PDF file
	model.PrintAllMarginalized(c_str(path + "out/Fit2nbbLV_plots.pdf"));

    // print all summary plots
	summary.PrintParameterPlot(c_str(path + "out/Fit2nbbLV_parameters.pdf"));
	summary.PrintCorrelationPlot(c_str(path + "out/Fit2nbbLV_correlation.pdf"));
	summary.PrintCorrelationMatrix(c_str(path + "out/Fit2nbbLV_correlationMatrix.pdf"));
    model.WriteHistosOnFile(path + "out/");
    // this will re-run the analysis without the LogLikelihood information
    BCLog::OutSummary("Building knowledge-update plots.");
    BCLog::SetLogLevelScreen(BCLog::warning);
	summary.PrintKnowledgeUpdatePlots(c_str(path + "out/Fit2nbbLV_update.pdf"));
    BCLog::SetLogLevelScreen(BCLog::summary);

    BCLog::OutSummary("Exiting");
	// close log file
	BCLog::CloseLog();

    if ( level != BCEngineMCMC::kLow ) {
        std::string command = "telegram-send \"runfit: completed\"";
        std::system(command.c_str());
    }

    return 0;
}
