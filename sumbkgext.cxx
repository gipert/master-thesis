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
    
    std::string phys;
    if      ( std::find(args.begin(), args.end(), "--K42" ) != args.end() ) phys = "K42";
    else if ( std::find(args.begin(), args.end(), "--K40" ) != args.end() ) phys = "K40";
    else    { std::cout << "Specify --K42 or --K40!\n"; return -1; }
// ----------------------------------------------------------------------------------------------------------
    
    // infos about runs
    GERDA::DataReader reader( std::string(std::getenv("GERDACPTDIR")) + "/misc/paths.txt", false, "GELATIO");
    // detector types
    GERDA::DetectorSet set("GELATIO");

    long int Ngen;
    double M;
    if ( phys == "K42" ) {
        Ngen = 5E09;
        M = 24563.385; // LAr mass [kg]
    }
    if ( phys == "K40" ) {
        Ngen = 1E08;
        M = 1.3615078; // fiber's mass [kg]
    }

    // fiber volume = 1.2966741074531858 E03 cm3
    // density 1.05 g/cm^3

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
    if ( phys == "K42" ) path = std::string(std::getenv("GERDACPTDIR")) + "/out/K42_homLAr_5E09.root";
    if ( phys == "K40" ) path = std::string(std::getenv("GERDACPTDIR")) + "/out/K40_onFiberShroud_1E08.root";
    TFile infile(path.c_str(), "READ"); 
    if (!infile.IsOpen()) { std::cerr << "Zombie infile!\n"; return -1; }

    for ( int i = 0; i < 40; ++i ) { 
        hist.push_back(dynamic_cast<TH1F*>(infile.Get(Form("h%i", i))));
        if (hist[i]->IsZombie()) { std::cout << "Zombie h" << std::to_string(i) << "!\n"; return -1; }
    }
    
    if ( phys == "K42" ) path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_homLAr.root";
    if ( phys == "K40" ) path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_K40onFiberShroud.root";
    TFile outfile(path.c_str(), "RECREATE");
    for ( int i = 0; i < 40; ++i ) {
        hist[i]->Scale(totalTime[i]*M/Ngen);
        hist[i]->Write();
    }

    TH1F histBEGe("energy_BEGe", "BEGe global MaGe energy spectrum", 7500, 0, 7500);
    TH1F histCOAX("energy_COAX", "COAX global MaGe energy spectrum", 7500, 0, 7500);
    TH1F histTotAll("energy_total", "global MaGe energy spectrum", 7500, 0, 7500);
       
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
