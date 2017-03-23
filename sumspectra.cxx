/* sumspectra.cxx
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

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <chrono>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"

#include "DataReader.h"
#include "ProgressBar.h"

int main( int argc, char** argv ) {

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
    std::vector<int> totalTime(40, 0);
    for ( const auto& i : dsm ) {
        for ( int j = 0; j < 40; ++j ) {
            if ( i.second[j] == 0 ) totalTime[j] += timeMap[i.first];
        }
    }

    // define reading objects
    std::unique_ptr<TFile> file;
    TTree* fTree;

    TTreeReader treereader;
    TTreeReaderArray<int>   det_id(treereader, "det_id");
    TTreeReaderArray<float> det_edep(treereader, "det_edep");

    // construct final histograms (MaGeOutput scheme because we are reading the MaGe output)
    std::vector<TH1F> hist;
    for ( int i = 0; i < 40; ++i ) {
        hist.emplace_back(Form("energy_det_tmp_id%i", i), Form("global MaGe energy spectrum, det_id = %i", i), 7500, 0, 7.5);
    }
    
    std::vector<TH1F> histTot;
    for ( int i = 0; i < 40; ++i ) {
        histTot.emplace_back(Form("energy_det_id%i", i), Form("global MaGe energy spectrum, det_id = %i", i), 7500, 0, 7.5);
    }


    ProgressBar bar;
    std::string filename;
    std::string display;
    int nentries;
    int size;
    float corrN;
    int corrTime;

// -----------------------------------------------------------------------------------------------------
    // lambda to fill histograms
    auto fillHistos = [&]( int i , std::string genopt , std::string phys ) {

        if      ( phys == "2nbbLV" ) filename = "/home/GERDA/pertoldi/simulations/2nbbLV/2nbbLV_";
        else if ( phys == "2nbb"   ) filename = "/home/GERDA/pertoldi/simulations/2nbb/2nbb_";
        display.clear();
        
        if ( genopt == "A_COAX" ) {
            filename += "AV_det11_";
            display += "AV_det11_";
            corrN = N76AV[i-1]/Ngen;
            corrTime = totalTime[i-1];
        }

        else if ( genopt == "D_COAX" ) {
            filename += "DV_det11_"; 
            display += "DV_det11_";
            corrN = N76DV[i-1]/Ngen;
            corrTime = totalTime[i-1];
        }

        else if ( genopt == "A_BEGe" ) {
            filename += "AV_det5_";
            display += "AV_det5_";
            corrN = N76AV[i+9]/Ngen;
            corrTime = totalTime[i+9];
        }
        
        else if ( genopt == "D_BEGe" ) {
            filename += "DV_det5_"; 
            display += "DV_det5_";
            corrN = N76DV[i+9]/Ngen;
            corrTime = totalTime[i+9];
        }

        filename += std::to_string(i) + ".root";
        display += std::to_string(i) + " ";
        
        file = std::unique_ptr<TFile>{ TFile::Open(filename.c_str(), "READ") };
        fTree = dynamic_cast<TTree*>(file->Get("fTree"));
        treereader.SetTree(fTree);       

        // fill with all entries
        if (verbose) { 
            nentries = treereader.GetEntries(true);
            bar.SetNIter(nentries);
            bar.Init(); 
        }
        int j = 0;
        std::cout << display << std::flush;
        while ( treereader.Next() ) {
            if (verbose) {bar.Update(j); j++;}
            size = det_id.GetSize();
            for ( int k = 0; k < size; ++k ) {
                if ( det_id[k] != 0 and det_id[k] != 1 and det_id[k] != 2 ) {
                    hist[det_id[k]].Fill(det_edep[k]);
                }
            }
        } std::cout << ' '; if (verbose) std::cout << nentries << " entries";
        
        // scale, add to final histogram, clear
        for ( int k = 0; k < 40; ++k ) {
            hist[i].Scale(corrN*corrTime);
            histTot[i].Add(&hist[i]);
            hist[i].Reset();
        }
        file->Close();
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
    
    for ( auto& h : hist ) {
        h.Write();
        histTotAll.Add(&h);
    }
    
    for ( int i = 0; i < 40; ++i ) {
        if ( i == 0 or i == 1 or i == 2 ) continue;
        else if ( i == 11 or i == 12 or i == 13 or i == 30 or i == 31 or i == 32 or i == 39 ) histCOAX.Add(&hist[i]);
        else histBEGe.Add(&hist[i]);
    }
    
    histCOAX.Write();
    histBEGe.Write();
    histTotAll.Write();

    fileout.Close();

    return 0;
}
