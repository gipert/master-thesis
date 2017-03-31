/* processbb.cxx
 *
 * Program to produce MaGe simulations for each detector (AV and DV)
 *
 * state-of-arts: enrCOAX are not summed, neither the events with
 * vertex into them nor energy deposition into them.
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 23/03/2017
 *
 */

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <chrono>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"

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
    // define reading objects
    std::unique_ptr<TFile> file;
    TTree* fTree;

    TTreeReader treereader;
    TTreeReaderArray<int>   det_id(treereader, "det_id");
    TTreeReaderArray<float> det_edep(treereader, "det_edep");

    // construct final histograms (MaGeOutput scheme because we are reading the MaGe output)
    std::vector<TH1F> hist;
    for ( int i = 0; i < 40; ++i ) {
        hist.emplace_back(Form("energy_det_id%i", i), Form("det_id = %i", i), 7500, 0, 7.5);
    }

    ProgressBar bar;
    std::string filename, outfilename;
    std::string display;
    int nentries;
    int size;

// -----------------------------------------------------------------------------------------------------
    // lambda to fill histograms
    auto fillHistos = [&]( int i , std::string genopt , std::string phys ) {

        if      ( phys == "2nbbLV" ) filename = "/home/GERDA/pertoldi/simulations/2nbbLV/2nbbLV_";
        else if ( phys == "2nbb"   ) filename = "/home/GERDA/pertoldi/simulations/2nbb/2nbb_";
        display.clear();

        if      ( phys == "2nbbLV" ) outfilename = "/home/GERDA/pertoldi/simulations/2nbbLV/processed/p_2nbbLV_";
        else if ( phys == "2nbb"   ) outfilename = "/home/GERDA/pertoldi/simulations/2nbb/processed/p_2nbb_";
        
        if ( genopt == "A_COAX" ) {
            filename += "AV_det11_";
            outfilename += "AV_det11_";
            display += "AV_det11_";
        }

        else if ( genopt == "D_COAX" ) {
            filename += "DV_det11_"; 
            outfilename += "DV_det11_";
            display += "DV_det11_";
        }

        else if ( genopt == "A_BEGe" ) {
            filename += "AV_det5_";
            outfilename += "AV_det5_";
            display += "AV_det5_";
        }
        
        else if ( genopt == "D_BEGe" ) {
            filename += "DV_det5_"; 
            outfilename += "DV_det5_";
            display += "DV_det5_";
        }

        filename += std::to_string(i) + ".root";
        outfilename += std::to_string(i) + ".root";
        display += std::to_string(i) + " ";
        
        file = std::unique_ptr<TFile>{ TFile::Open(filename.c_str(), "READ") };
        fTree = dynamic_cast<TTree*>(file->Get("fTree"));
        treereader.SetTree(fTree);       

        // fill with all entries
        if (verbose) { 
            nentries = treereader.GetEntries(true);
            bar.SetNIter(nentries);
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
        
        // let's save on disk
        TFile outfile(outfilename.c_str(), "RECREATE");
        for ( auto& h : hist ) { 
            h.Write();
            h.Reset();
        }
        outfile.Close();
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

    return 0;
}
