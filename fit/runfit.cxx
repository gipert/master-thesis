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

double GetPValue(Fit2nbbLV& m, BCEngineMCMC::Precision level, bool save = false);

int main( int argc, char** argv ) {

/////////////////////////////////////////////
    const int rangeUp = 5300;  // [keV]
    const int rangeDown = 570; // [keV] above 39Ar Q-value
    BCEngineMCMC::Precision level(BCEngineMCMC::kLow);
/////////////////////////////////////////////

    auto c_str = [](std::string s) { return s.c_str(); };

    // save command line arguments
    std::vector<std::string> args;
    for ( int i = 0; i < argc; ++i ) args.emplace_back(argv[i]);

    // help
    if ( std::find(args.begin(), args.end(), "--help") != args.end() ) {
        std::cout << std::endl
                  << "Usage:\n\n"
                  << "    runfit [OPTIONS]\n\n"
                  << "Options:\n\n"
                  << "    --fixbinning     : use fixed-size binning instead of the\n"
                  << "                       default one (variable)\n"
                  << "    --outdir [DIR]   : set directory to store results\n"
                  << "                       default : $GERDACPTDIR/out\n"
                  << "    --fixfile [FILE] : set fixfile\n\n";
        return 0;
    }

    // retrieve the name of the directory containing the output
    auto result = std::find( args.begin(), args.end(), "--outdir");
    std::string outdirname;
    if ( result != args.end() ) outdirname = *(result+1);
    else                        outdirname = std::string(std::getenv("GERDACPTDIR")) + "/out";
    std::system(c_str("mkdir -p " + outdirname));

    // retrieve the name of the fix file
    result = std::find( args.begin(), args.end(), "--fixfile");
    std::string fixfilepath;
    if ( result != args.end() ) fixfilepath = *(result+1);
    else                        fixfilepath = std::string(std::getenv("GERDACPTDIR")) + "/misc/fallback.fix";

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
    // [14] Bi212holder
    // [15] Tl208holder
    // [16] Pb214holder
    // [17] Bi214holder
    // [18] K40cables
    // [19] Bi212cables
    // [20] Tl208cables
    // [21] Pb214cables
    // [22] Bi214cables
    // [23] K40minishroud
    // [24] Pa234minishroud
    // [25] Bi207minishroud
    // [26] Bi207cables
    // [27] Bi207holder
    // [28] Pb214minishroud
    // [29] Bi214minishroud
    // [30] K42minishroudsurface
    // [31] Pa234cables
    // [32] Pa234holders
    // [33] K42homLArAA

    std::vector<std::string> simpath;

    /*[0]*/  simpath.push_back(path + "sumMaGe_2nbb.root");
    /*[1]*/  simpath.push_back(path + "sumMaGe_2nbbLV.root");
    /*[2]*/  simpath.push_back(path + "sumMaGe_K42homLAr.root");
    /*[3]*/  simpath.push_back(path + "sumMaGe_K40fibers.root");
    /*[4]*/  simpath.push_back(path + "sumMaGe_Bi212fibers.root");
    /*[5]*/  simpath.push_back(path + "sumMaGe_Tl208fibers.root");
    /*[6]*/  simpath.push_back(path + "sumMaGe_Pb214fibers.root");
    /*[7]*/  simpath.push_back(path + "sumMaGe_Bi214fibers.root");
    /*[8]*/  simpath.push_back(path + "alpha_model_run53-79.root");
    /*[9]*/  simpath.push_back(path + "sumMaGe_K42nPlus.root");
    /*[10]*/ simpath.push_back(path + "sumMaGe_K42pPlus.root");
    /*[11]*/ simpath.push_back(path + "sumMaGe_Ac228holder.root");
    /*[12]*/ simpath.push_back(path + "sumMaGe_Co60holder.root");
    /*[13]*/ simpath.push_back(path + "sumMaGe_K40holder.root");
    /*[14]*/ simpath.push_back(path + "sumMaGe_Bi212holder.root");
    /*[15]*/ simpath.push_back(path + "sumMaGe_Tl208holder.root");
    /*[16]*/ simpath.push_back(path + "sumMaGe_Pb214holder.root");
    /*[17]*/ simpath.push_back(path + "sumMaGe_Bi214holder.root");
    /*[18]*/ simpath.push_back(path + "sumMaGe_K40cables.root");
    /*[19]*/ simpath.push_back(path + "sumMaGe_Bi212cables.root");
    /*[20]*/ simpath.push_back(path + "sumMaGe_Tl208cables.root");
    /*[21]*/ simpath.push_back(path + "sumMaGe_Pb214cables.root");
    /*[22]*/ simpath.push_back(path + "sumMaGe_Bi214cables.root");
    /*[23]*/ simpath.push_back(path + "sumMaGe_K40minishroud.root");
    /*[24]*/ simpath.push_back(path + "sumMaGe_Pa234minishroud.root");
    /*[25]*/ simpath.push_back(path + "sumMaGe_Bi207minishroud.root");
    /*[26]*/ simpath.push_back(path + "sumMaGe_Bi207cables.root");
    /*[27]*/ simpath.push_back(path + "sumMaGe_Bi207holder.root");
    /*[28]*/ simpath.push_back(path + "sumMaGe_Pb214minishroud.root");
    /*[29]*/ simpath.push_back(path + "sumMaGe_Bi214minishroud.root");
    /*[30]*/ simpath.push_back(path + "sumMaGe_K42minishroudsurface.root");
    /*[31]*/ simpath.push_back(path + "sumMaGe_Pa234cables.root");
    /*[32]*/ simpath.push_back(path + "sumMaGe_Pa234holder.root");
    /*[33]*/ simpath.push_back(path + "sumMaGe_K42homLArAA.root");

    // TODO: push_back new files here

    for ( const auto& p : simpath ) simFile.emplace_back( new TFile(p.c_str(), "READ") );

    for ( auto& f : simFile ) if (!f->IsOpen()) { std::cout << "At least one zombie simFile!\n"; return -1; }

    TH1D* hDataBEGetmp;   fileData.GetObject("energyBEGeAll",    hDataBEGetmp);
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

    // remove events in blinded window Qbb +- 25 keV
    for ( int i = 2014; i <= 2064; i++ ) {
        hDataBEGetmp->SetBinContent(i,0);
        hDataCOAXtmp->SetBinContent(i,0);
        for ( int j = 0; j < (int)hSimBEGetmp.size(); ++j ) {
            hSimBEGetmp[j]->SetBinContent(i,0);
            hSimCOAXtmp[j]->SetBinContent(i,0);
        }
    }

    // rebin the histograms in new histograms
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
    std::vector<int>                 vDataCOAX(nBins),          vDataBEGe(nBins);
    std::vector<std::vector<double>> vSimCOAX(hSimCOAX.size()), vSimBEGe(hSimBEGe.size());

    // alphas need to be normalized, we want to measure counts
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
    path = outdirname + "/";
    BCModelOutput output(&model, c_str(path + "markowChains.root"));
    model.WriteMarkovChain(true);

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::SetLogLevelFile(BCLog::detail);
    BCLog::SetLogLevelScreen(BCLog::summary);
    BCLog::OpenLog(c_str(path + "logBAT.txt"));

    // set precision (number of samples in Markov chain)
    model.MCMCSetPrecision(level);

    // set data
    model.SetBinning(dbin);
    model.SetDataBEGe(vDataBEGe);
    model.SetDataCOAX(vDataCOAX);
    model.SetSimBEGe(vSimBEGe);
    model.SetSimCOAX(vSimCOAX);
    model.SetFitRange(rangeDown, rangeUp);

    BCLog::OutSummary(c_str("Saving results in " + path));
    if ( std::find( args.begin(), args.end(), "--fixbinning" ) == args.end() ) BCLog::OutSummary("Adopting variable binning size ");
    else BCLog::OutSummary("Adopting fixed binning size");
    // eventually fix parameters as indicated in external file
    BCLog::OutSummary(c_str("Fix file: " + fixfilepath));
    std::ifstream fixfile(fixfilepath);
    int n, p;
    std::string comment;
    while ( fixfile >> p >> n >> comment ) {
        if ( n == 0 ) {
            model.GetParameter(p)->Fix(0);
            std::cout << "Summary : Fixing parameter: " << model.GetParameter(p)->GetName() << '\n';
        }
    }
    // set parameter binning
    //m.SetNbins(1000);

    // run MCMC and marginalize posterior w/r/t all parameters and all
    // combinations of two parameters
    auto start = std::chrono::system_clock::now();
    model.MarginalizeAll();
    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now()-start);
    std::cout << "Summary : Time spent: " << elapsed.count() << "s.\n";

    // run mode finding, by default using Minuit, working shitty
    //BCLog::SetLogLevelScreen(BCLog::detail);
    //model.FindMode(model.GetBestFitParameters());
    //BCLog::SetLogLevelScreen(BCLog::summary);

    std::cout << std::endl;
    double pvalue = GetPValue(model, level, false);
    std::cout << "Summary : pValue = " << pvalue << std::endl;

    // OUTPUT
    // print results of the analysis into a text file
    model.PrintResults(c_str(path + "Fit2nbbLV_results.txt"));
    // draw all marginalized distributions into a PDF file
    //model.PrintAllMarginalized(c_str(path + "Fit2nbbLV_plots.pdf"));

    // print all summary plots
    summary.PrintParameterPlot(c_str(path + "Fit2nbbLV_parameters.pdf"));
    summary.PrintCorrelationPlot(c_str(path + "Fit2nbbLV_correlation.pdf"));
    summary.PrintCorrelationMatrix(c_str(path + "Fit2nbbLV_correlationMatrix.pdf"));
    model.WriteHistosOnFile(path);
    // this will re-run the analysis without the LogLikelihood information
    BCLog::OutSummary("Building knowledge-update plots.");
    BCLog::SetLogLevelScreen(BCLog::warning);
    summary.PrintKnowledgeUpdatePlots(c_str(path + "Fit2nbbLV_update.pdf"));
    BCLog::SetLogLevelScreen(BCLog::summary);

    BCLog::OutSummary("Exiting");
    // close log file
    BCLog::CloseLog();

    return 0;
}
