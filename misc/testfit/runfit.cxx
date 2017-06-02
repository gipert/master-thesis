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
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <memory>

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCModelOutput.h>
#include <BAT/BCParameter.h>

#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"

int main( int argc, char** argv ) {

/////////////////////////////////////////////
    const int rangeUp = 2039;  // [keV]
    const int rangeDown = 570; // [keV] above 39Ar Q-value
    BCEngineMCMC::Precision level(BCEngineMCMC::kMedium);
/////////////////////////////////////////////

    auto c_str = [](std::string s) { return s.c_str(); };

    // save command line arguments
    std::vector<std::string> args;
    for ( int i = 0; i < argc; ++i ) args.emplace_back(argv[i]);

    TH1::AddDirectory(false);
    // retrieve data
    std::string path = std::string(std::getenv("GERDACPTDIR")) + "/data/";
    TFile fileData(c_str("../sim2nbbspectrum.root"), "READ");
    if (!fileData.IsOpen()) { std::cout << "Zombie fileData!\n"; return -1; }

    TH1F* hDataBEGetmp; fileData.GetObject("energyBEGe", hDataBEGetmp);
    if (!hDataBEGetmp) std::cout << "zombie\n";

    std::vector<TH1*> hSimBEGetmp;

    TH1* tmp;

    TFile * f = new TFile(c_str(path + "sumMaGe_2nbb.root"), "READ");
    f->GetObject("energy_BEGe", tmp);
    hSimBEGetmp.push_back(tmp);

    f = new TFile(c_str( path + "sumMaGe_2nbbLV.root"), "READ");
    f->GetObject("energy_BEGe", tmp);
    hSimBEGetmp.push_back(tmp);

    fileData.Close();

    for ( auto& h : hSimBEGetmp ) if (!h) { std::cout << "There's at least one zombie simBEGe hist!\n"; return -1; }
    // create binning

    int nBins;
    std::vector<double> dbin; // this is what ROOT wants

    // variable
    if ( std::find( args.begin(), args.end(), "--fixbinning" ) == args.end() ) {
        nBins = 1847;
        dbin = std::vector<double>(nBins+1);
        std::vector<int> avoid = { 568, 572, 580, 584, 608, 612, 908,
                                   912, 968, 972, 1000, 1004, 1060, 1064,
                                   1120, 1172, 1176, 1236, 1240, 1332, 1460,
                                   1464, 1524, 1528, 1764, 2204, 2612, 2616 };

        int k = 0, i = 0;
        while (1) {
            if ( std::find( avoid.begin(), avoid.end(), k) == avoid.end() ) {
                dbin[i] = k;
                i++;
            }
            k += 4;
            if ( k > 7500 ) break;
        }
    }

    // fixed
    else {
        nBins = 1875;
        dbin = std::vector<double>(nBins+1);
        int k = 0;
        for ( int i = 0; i <= nBins; ++i ) {
            dbin[i] = k;
            k+=4;
        }
    }

    // rebin the histograms in new histograms
    TH1F* hDataBEGe = dynamic_cast<TH1F*>(hDataBEGetmp->Rebin(nBins, "hDataBEGetmp", &dbin[0]));

    std::vector<TH1*> hSimBEGe;

    for ( unsigned int i = 0; i < hSimBEGetmp.size(); ++i ) {
       hSimBEGe.push_back(dynamic_cast<TH1F*>(hSimBEGetmp[i]->Rebin(nBins, Form("hSimBEGe%i", i), &dbin[0])));
    }

    delete hDataBEGetmp;
    hSimBEGetmp.clear();

    // convert in std::vector
    std::vector<int>                 vDataCOAX(nBins),          vDataBEGe(nBins);
    std::vector<std::vector<double>> vSimBEGe(hSimBEGe.size());


    for ( int i = 0; i < nBins; ++i ) {
        vDataBEGe[i] = hDataBEGe->GetBinContent(i+1);

        for ( unsigned int j = 0; j < hSimBEGe.size(); ++j ) {
            vSimBEGe[j].push_back(hSimBEGe[j]->GetBinContent(i+1));
        }
    }

    delete hDataBEGe;
    hSimBEGe.clear();

// ==========================================================================================================

    // create Fit2nbbLV object
    Fit2nbbLV model("Fit2nbbLV");
    // create a new summary tool object
    BCSummaryTool summary(&model);

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::SetLogLevelFile(BCLog::detail);
    BCLog::SetLogLevelScreen(BCLog::detail);
    BCLog::OpenLog(c_str("logBAT.txt"));

    // set precision (number of samples in Markov chain)
    model.MCMCSetPrecision(level);

    // set data
    model.SetBinning(dbin);
    model.SetDataBEGe(vDataBEGe);
    model.SetSimBEGe(vSimBEGe);
    model.SetFitRange(rangeDown, rangeUp);

    if ( std::find( args.begin(), args.end(), "--fixbinning" ) == args.end() ) BCLog::OutSummary("Adopting variable binning size ");
    else BCLog::OutSummary("Adopting fixed binning size");

    // run MCMC and marginalize posterior w/r/t all parameters and all
    // combinations of two parameters
    model.SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
    model.MarginalizeAll();

    // OUTPUT
    // print results of the analysis into a text file
    model.PrintResults(c_str("Fit2nbbLV_results.txt"));
    // draw all marginalized distributions into a PDF file
    //model.PrintAllMarginalized(c_str(path + "Fit2nbbLV_plots.pdf"));

    // print all summary plots
    summary.PrintParameterPlot(c_str("Fit2nbbLV_parameters.pdf"));
    summary.PrintCorrelationPlot(c_str("Fit2nbbLV_correlation.pdf"));
    summary.PrintCorrelationMatrix(c_str("Fit2nbbLV_correlationMatrix.pdf"));
    // this will re-run the analysis without the LogLikelihood information
    BCLog::OutSummary("Building knowledge-update plots.");
    BCLog::SetLogLevelScreen(BCLog::warning);
    summary.PrintKnowledgeUpdatePlots(c_str("Fit2nbbLV_update.pdf"));
    BCLog::SetLogLevelScreen(BCLog::summary);

    BCLog::OutSummary("Exiting");
    // close log file
    BCLog::CloseLog();

    return 0;
}
