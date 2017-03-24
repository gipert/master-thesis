/* sumMaGe.cxx
 *
 * Program to sum MaGe simulations for each detector (AV and DV)
 * into global energy spectra.
 *
 * state-of-arts: enrCOAX are not summed, neither the events with
 * vertex into them nor energy deposition into them. Different
 * volumes and live times of the detectors (depending on the 
 * considered runs) are taken into account.
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
#include <memory>

#include "TFile.h"
#include "TH1F.h"

#include "DataReader.h"

int main( int argc, char** argv ) {
    
    TH1::AddDirectory(kFALSE);

    std::vector<std::string> args(argc);
    for ( int i = 0; i < argc; ++i ) args[i] = argv[i];
    
    bool verbose = false;
    if ( std::find(args.begin(), args.end(), "--verbose"  ) != args.end() ) verbose = true;
    
    std::string phys;
    if      ( std::find(args.begin(), args.end(), "--2nbb"  ) != args.end() ) phys = "2nbb";
    else if ( std::find(args.begin(), args.end(), "--2nbbLV") != args.end() ) phys = "2nbbLV";
    else { std::cout << "Please specify --2nbb or --2nbbLV option!\n"; return 0; }
// ----------------------------------------------------------------------------------------------------------
    // infos about runs
    GERDA::DataReader reader( std::string(std::getenv("GERDACPTDIR")) + "/misc/paths.txt", false, "MaGeInput");
    // infos on experimental setup
    GERDA::DetectorSet set("MaGeInput");
    
    ////// get N76 with MaGe input naming convention [mol]
    std::vector<float> N76AV = set.GetActiveN76Ge();
    std::vector<float> N76DV = set.GetDeadN76Ge(); 
    // number of generated events
    int Ngen = 10E07;

    // get detectorStatusMap for selected runs
    auto dsm = reader.GetDetectorStatusMap();

    // get live times saved in out/results.dat and store [yr]
    std::string path = std::string(std::getenv("GERDACPTDIR")) + "/out/results.dat";
    std::ifstream timeFile(path.c_str());
    std::map<int,double> timeMap;
    unsigned int runID, time;
    while ( timeFile >> runID >> time ) timeMap.insert(std::make_pair(runID,time*1./31536000));

    ////// calculate the total time each detector is ON --> detector status = 0
    std::vector<double> totalTime(40, 0);
    for ( const auto& i : dsm ) {
        for ( int j = 0; j < 40; ++j ) {
            if ( i.second[j] == 0 ) totalTime[j] += timeMap[i.first];
        }
    }
    
    // construct histograms (MaGeOutput scheme because we are reading the MaGe output)
    std::vector<std::unique_ptr<TH1F>> hist;
    // final histogram
    std::vector<TH1F> histTot;
    for ( int i = 0; i < 40; ++i ) {
        histTot.emplace_back(Form("energytot_det_id%i", i), Form("det_id = %i", i), 7500, 0, 7.5);
    }

    std::string filename;
    double corrN = 0, corrTime = 0;

// -----------------------------------------------------------------------------------------------------
    // lambda to fill histograms
    auto fillHistos = [&]( int i , std::string genopt , std::string phys ) {

        if      ( phys == "2nbbLV" ) filename = "/home/GERDA/pertoldi/simulations/2nbbLV/processed/p_2nbbLV_";
        else if ( phys == "2nbb"   ) filename = "/home/GERDA/pertoldi/simulations/2nbb/processed/p_2nbb_";
        
        if ( genopt == "A_COAX" ) {
            filename += "AV_det11_";
            corrN = N76AV[i-1]/Ngen;
            corrTime = totalTime[i-1];
        }

        else if ( genopt == "D_COAX" ) {
            filename += "DV_det11_"; 
            corrN = N76DV[i-1]/Ngen;
            corrTime = totalTime[i-1];
        }

        else if ( genopt == "A_BEGe" ) {
            filename += "AV_det5_";
            corrN = N76AV[i+9]/Ngen;
            corrTime = totalTime[i+9];
        }
        
        else if ( genopt == "D_BEGe" ) {
            filename += "DV_det5_"; 
            corrN = N76DV[i+9]/Ngen;
            corrTime = totalTime[i+9];
        }

        filename += std::to_string(i) + ".root";
        
        // open file
        TFile file(filename.c_str(), "READ");
        
        // retrieve histograms
        for ( int j = 0; j < 40; ++j ) {
            hist.emplace_back( dynamic_cast<TH1F*>(file.Get(Form("energy_det_id%i", j))) );
        }
        // scale, add to final histogram, clear
        for ( int k = 0; k < 40; ++k ) {
            hist[k]->Scale(corrN*corrTime);
            histTot[k].Add(hist[k].get());
        }
        hist.clear();
        file.Close();
        return;
    };
// -----------------------------------------------------------------------------------------------------
    
    // loop over enrCOAX files
    for ( int i = 4; i <= 10; ++i ) {
        
        auto start = std::chrono::system_clock::now();
        
        // run
        fillHistos(i, "A_COAX", phys); if (verbose) std::cout << std::endl;
        fillHistos(i, "D_COAX", phys);
        
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
        if (verbose) std::cout << " [" << elapsed.count()*1./1000 << "s]\n";
    }

    // loop over BEGe
    for ( int i = 1; i <= 30; ++i ) {

        auto start = std::chrono::system_clock::now();
        
        // run
        fillHistos(i, "A_BEGe", phys); if (verbose) std::cout << std::endl;
        fillHistos(i, "D_BEGe", phys);

        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
        if (verbose) std::cout << " [" << elapsed.count()*1./1000 << "s]\n";
    }

    TH1F histBEGe("energy_BEGe", "BEGe global MaGe energy spectrum", 7500, 0, 7.5);
    TH1F histCOAX("energy_COAX", "COAX global MaGe energy spectrum", 7500, 0, 7.5);
    TH1F histTotAll( "energy_total", "global MaGe energy spectrum", 7500, 0, 7.5 );
    
    if ( phys == "2nbb" ) path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_2nbb.root";
    else path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumMaGe_2nbbLV.root";
    TFile fileout(path.c_str(), "RECREATE");
    
    for ( auto& h : histTot ) {
        h.Write();
        histTotAll.Add(&h);
    }
    
    for ( int i = 0; i < 40; ++i ) {
        if ( i == 0 or i == 1 or i == 2 ) continue;
        else if ( i == 11 or i == 12 or i == 13 or 
                  i == 30 or i == 31 or i == 32 or i == 39 ) histCOAX.Add(&histTot[i]);
        else histBEGe.Add(&histTot[i]);
    }
    
    histCOAX.Write();
    histBEGe.Write();
    histTotAll.Write();

    fileout.Close();

    return 0;
}
