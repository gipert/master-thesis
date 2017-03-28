/* sumbbkgext.cxx
 *
 * Program to sum MaGe simulations for external bkg sources
 * into global energy spectra.
 *
 * supported sources:
 *   homLAr
 *
 * state-of-arts: enrCOAX are not summed. Different live times 
 * of the detectors (depending on the considered runs) are taken 
 * into account.
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 15/03/2017
 *
 */

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <map>

#include "TFile.h"
#include "TH1F.h"

#include "DataReader.h"

int main( int argc, char** argv ) {
    
    std::vector<std::string> args(argc);
    for ( int i = 0; i < argc; ++i ) args[i] = argv[i];
    
    //bool verbose = false;
    //if ( std::find(args.begin(), args.end(), "--verbose"  ) != args.end() ) verbose = true; 
// ----------------------------------------------------------------------------------------------------------
    
    // infos about runs
    GERDA::DataReader reader( std::string(std::getenv("GERDACPTDIR")) + "/misc/paths.txt", false, "GELATIO");
    // detector types
    GERDA::DetectorSet set("GELATIO");

    long int Ngen = 5E09;
    // LAr mass [kg]
    double M = 24563.385;

    // get detectorStatusMap for selected runs
    auto dsm = reader.GetDetectorStatusMap();

    // get live times saved in out/results.dat and store [s]
    std::string path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumData.dat";
    std::ifstream timeFile(path.c_str());
    std::map<int,int> timeMap;
    int runID, time;
    while ( timeFile >> runID >> time ) timeMap.insert(std::make_pair(runID,time));

    ////// calculate the total time each detector is ON --> detector status = 0
    std::vector<int> totalTime(40, 0);
    for ( const auto& i : dsm ) {
        for ( int j = 0; j < 40; ++j ) {
            if ( i.second[j] == 0 ) totalTime[j] += timeMap[i.first];
        }
    }
    
    std::vector<TH1F*> hist;

    // get original spectra
    path = std::string(std::getenv("GERDACPTDIR")) + "/out/K42_homLAr_5E09.root";
    TFile infile(path.c_str(), "READ"); 
    if (!infile.IsOpen()) { std::cerr << "Zombie infile!\n"; return -1; }

    for ( int i = 0; i < 40; ++i ) { 
        hist.push_back(dynamic_cast<TH1F*>(infile.Get(Form("h%i", i))));
        if (hist[i]->IsZombie()) { std::cout << "Zombie h" << std::to_string(i) << std::endl; return -1; }
    }
    
    path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_homLAr.root";
    TFile outfile(path.c_str(), "RECREATE");
    for ( int i = 0; i < 40; ++i ) {
        hist[i]->Scale(totalTime[i]*M/Ngen);
        hist[i]->Write();
    }

    TH1F histBEGe("energy_BEGe", "BEGe global MaGe energy spectrum", 8500, 0, 8500);
    TH1F histCOAX("energy_COAX", "COAX global MaGe energy spectrum", 8500, 0, 8500);
    TH1F histTotAll("energy_total", "global MaGe energy spectrum", 8500, 0, 8500);
       
    for ( auto& h : hist ) histTotAll.Add(h);
    
    for ( int i = 0; i < 40; ++i ) {
        if      ( set.GetDetectorTypes()[i] == 3 ) continue;
        else if ( set.GetDetectorTypes()[i] == 2 ) histCOAX.Add(hist[i]);
        else                                       histBEGe.Add(hist[i]);
    }
    
    histCOAX.Write();
    histBEGe.Write();
    histTotAll.Write();

    outfile.Close();
    infile.Close();

    return 0;
}
