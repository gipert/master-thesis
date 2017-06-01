/* sim2nbbspectrum.cxx
 *
 * Program to simulate energy spectra (BEGe and EnrCoax)
 *
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 01/05/2017
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
#include "DetectorSet.h"

int main( int argc, char** argv ) {

    std::vector<std::string> args(argc);
    for ( int i = 0; i < argc; ++i ) args[i] = argv[i];

    bool verbose = false;
    if ( std::find(args.begin(), args.end(), "--verbose"  ) != args.end() ) verbose = true;

    std::string phys;
    if      ( std::find(args.begin(), args.end(), "--2nbb"  ) != args.end() ) phys = "2nbb";
    else if ( std::find(args.begin(), args.end(), "--2nbbLV") != args.end() ) phys = "2nbbLV";
    else { std::cout << "Please specify --2nbb or --2nbbLV option!\n"; return 0; }

    // retrieve number of 2nbb events
    int N2nbb;
    auto result = std::find(args.begin(), args.end(), "--Nev" );
    if ( result != args.end() ) N2nbb = std::stoi(*(result+1));
    else { std::cout << "Please specify the number of 2nbb events!\n"; return -1; }
// ----------------------------------------------------------------------------------------------------------

    // define reading objects
    std::unique_ptr<TFile> file;
    TTree* fTree;

    TTreeReader treereader;
    TTreeReaderArray<int>   det_id(treereader, "det_id");
    TTreeReaderArray<float> det_edep(treereader, "det_edep");

    // detectors' properties
    GERDA::DetectorSet set("MaGeInput");

    TH1F histBEGe("energyBEGe", "energyBEGe", 7500, 0, 7500);
    //TH1F histCOAX("energyCOAX", "energyCOAX", 7500, 0, 7500);

    ProgressBar bar;
    std::string filename;
    std::string display;
    int nentries;
    int size;

    // determine number of events from each volume
    std::vector<int> AV_Nev;
    std::vector<int> DV_Nev;

    // compute total volume (BEGe);
    double totVol = 0;
    for ( int i = 0; i < 40; ++i ) if ( set.GetDetectorTypes()[i] == 1 ) totVol += set.GetVolume()[i];

    for ( int i = 0; i < 40; ++i ) {
        if ( set.GetDetectorTypes()[i] == 1 ) {
            AV_Nev.push_back(((double)set.GetActiveVolume()[i]/totVol)*N2nbb);
            DV_Nev.push_back(((double)set.GetDeadVolume()[i]/totVol)*N2nbb);
        }
    }

// -----------------------------------------------------------------------------------------------------
    // lambda to fill histograms
    auto fillHistos = [&]( int i , std::string genopt , std::string phys ) {

        std::string basename = "/storage/gpfs_data/gerda/gerda_scratch/pertoldi/simulations/";

        if      ( phys == "2nbbLV" ) filename = basename + "2nbbLV/2nbbLV_";
        else if ( phys == "2nbb"   ) filename = basename + "2nbb/2nbb_";
        display.clear();

        if ( genopt == "A_COAX" ) {
            filename += "AV_det11_";
            display += "AV_det11_";
        }

        else if ( genopt == "D_COAX" ) {
            filename += "DV_det11_"; 
            display += "DV_det11_";
        }

        else if ( genopt == "A_BEGe" ) {
            filename += "AV_det5_";
            display += "AV_det5_";
            nentries = AV_Nev[i-1];
        }

        else if ( genopt == "D_BEGe" ) {
            filename += "DV_det5_"; 
            display += "DV_det5_";
            nentries = DV_Nev[i-1];
        }

        filename += std::to_string(i) + ".root";
        display += std::to_string(i) + ": " + std::to_string(nentries) + " events: ";

        file = std::unique_ptr<TFile>{ TFile::Open(filename.c_str(), "READ") };
        fTree = dynamic_cast<TTree*>(file->Get("fTree"));
        treereader.SetTree(fTree);

        // fill with all entries
        if (verbose) {
            bar.SetNIter(nentries);
        }
        int j = 0;
        std::cout << display << std::flush;
        for ( int i = 0; i < nentries; ++i ) {
            treereader.Next();
            if (verbose) {bar.Update(j); j++;}
            size = det_id.GetSize();
            for ( int k = 0; k < size; ++k ) {
                // NOTE: ignoring energy depositions in GTFs and GD02D
                if ( det_id[k] != 0 and det_id[k] != 1 and det_id[k] != 2 and det_id[k] != 9 ) {
                    if ( genopt == "A_BEGe" or genopt == "D_BEGe" ) histBEGe.Fill(det_edep[k]*1000);
                    //if ( genopt == "A_COAX" or genopt == "D_COAX" ) histCOAX.Fill(det_edep[k]*1000);
                }
            }
        } std::cout << ' '; if (verbose) std::cout << nentries << " entries";

        file->Close();
        std::cout << "done. \n";
        return;
    };
// -----------------------------------------------------------------------------------------------------
/*
    // loop over COAX files
    for ( int i = 1; i <= 10; ++i ) {

        auto start = std::chrono::system_clock::now();

        // run
        fillHistos(i, "A_COAX", phys); if (verbose) std::cout << std::endl;
        fillHistos(i, "D_COAX", phys);

        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
        if (verbose) std::cout << " [" << elapsed.count()*1./1000 << "s]\n";
    }
*/
    // loop over BEGe
    for ( int i = 1; i <= 30; ++i ) {

        auto start = std::chrono::system_clock::now();

        // run
        fillHistos(i, "A_BEGe", phys); if (verbose) std::cout << std::endl;
        fillHistos(i, "D_BEGe", phys);

        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
        if (verbose) std::cout << " [" << elapsed.count()*1./1000 << "s]\n";
    }

    // save on disk
    std::string outname(std::getenv("GERDACPTDIR"));
    if ( phys == "2nbb" ) outname += "/misc/sim2nbbspectrum.root";
    if ( phys == "2nbbLV" ) outname += "/misc/sim2nbbLVspectrum.root";
    TFile outfile( outname.c_str() ,"RECREATE");

    histBEGe.Write();
    //histCOAX.Write();

    return 0;
}
