/* 
 * Program to sum MaGe simulations for each detector (AV and DV)
 * into one global energy spectra.
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 15/03/2017
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

#include "DataReader.h"
#include "ProgressBar.h"

int main() {

    // get configs
    std::ifstream input("misc/paths.txt");
    if ( !input.is_open() ) { std::cerr << "File with paths not found! Aborting...\n"; return 0; }
    std::string metapath, datapath, configpath;
    input >> metapath >> datapath >> configpath;
    input.close();

    GERDA::DataReader reader(metapath, datapath, configpath);
    
    // get volumes
    std::vector<float> AV = reader.GetActiveVolume("MaGe");
    std::vector<float> DV = reader.GetDeadVolume("MaGe");
  
    // find max volume
    auto maxvolumeAV = std::max_element(AV.begin()+3, AV.end());
    auto maxvolumeDV = std::max_element(DV.begin()+3, DV.end());
    float maxvolume = *maxvolumeAV > *maxvolumeDV ? *maxvolumeAV : *maxvolumeDV;

    // define reading objects
    std::unique_ptr<TFile> file;
    TTree* fTree;

    TTreeReader treereader;
    TTreeReaderArray<int>   det_id(treereader, "det_id");
    TTreeReaderArray<float> det_edep(treereader, "det_edep");

    // final histogram
    std::vector<TH1F> hist;
    for ( int i = 0; i < 40; i++ ) {
        hist.emplace_back(Form("energy_det_id%i", i), Form("global MaGe energy spectrum, det_id = %i", i), 2100, 0, 2.1);
    }
    
    ProgressBar bar;
    std::string filename;
    std::string display;
    int nentries;
    int size;
    float corrFactor;

// -----------------------------------------------------------------------------------------------------
    // lambda to fill histograms
    auto fillHistos = [&]( int i , std::string opt ) {
        
        if ( opt == "A_COAX" ) {
            filename = "/home/GERDA/pertoldi/simulations/2nbbLV/2nbbLV_AV_det11_"  
                       + std::to_string(i) + ".root";
            display = "AV_det11_" + std::to_string(i) + " ";
            corrFactor = AV[i-1]/maxvolume;
        }

        else if ( opt == "D_COAX" ) {
            filename = "/home/GERDA/pertoldi/simulations/2nbbLV/2nbbLV_DV_det11_" 
                       + std::to_string(i) + ".root";
            display = "DV_det11_" + std::to_string(i) + " ";
            corrFactor = DV[i-1]/maxvolume;
        }

        else if ( opt == "A_BEGe" ) {
            filename = "/home/GERDA/pertoldi/simulations/2nbbLV/2nbbLV_AV_det5_" 
                       + std::to_string(i) + ".root";
            display = "AV_det5_" + std::to_string(i) + " ";
            corrFactor = AV[i+9]/maxvolume;
        }
        
        else if ( opt == "D_BEGe" ) {
            filename = "/home/GERDA/pertoldi/simulations/2nbbLV/2nbbLV_DV_det5_" 
                       + std::to_string(i) + ".root";
            display = "DV_det5_" + std::to_string(i) + " ";
            corrFactor = DV[i+9]/maxvolume;
        }

        else { std::cout << "wut?\n"; return; }
        
        file = std::unique_ptr<TFile>{ TFile::Open(filename.c_str(), "READ") };
        fTree = dynamic_cast<TTree*>(file->Get("fTree"));
        treereader.SetTree(fTree);       

        // fill with correct entries
        nentries = treereader.GetEntries(true)*corrFactor;
        bar.SetNIter(nentries);
        std::cout << display;
        bar.Init();
        for ( int j = 0; j < nentries; j++ ) {
            bar.Update(j);
            treereader.Next();
            size = det_id.GetSize();
            for ( int k = 0; k < size; k++ ) {
                if ( det_id[k] != 0 and det_id[k] != 1 and det_id[k] != 2 ) {
                    hist[det_id[k]].Fill(det_edep[k]);
                }
            }
        } std::cout << ' ' << nentries << " entries";
        file->Close();
        return;
    };
// -----------------------------------------------------------------------------------------------------
    
    // loop over enrCOAX files
    for ( int i = 4; i <= 10; i++ ) {
        
        auto start = std::chrono::system_clock::now();

        fillHistos(i, "A_COAX"); std::cout << std::endl;
        fillHistos(i, "D_COAX");
        
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
        std::cout << " [" << elapsed.count()*1./1000 << "s]\n";
    }

    // loop over BEGe
    for ( int i = 1; i <= 30; i++ ) {

        auto start = std::chrono::system_clock::now();
        
        fillHistos(i, "A_BEGe"); std::cout << std::endl;
        fillHistos(i, "D_BEGe");

        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
        std::cout << " [" << elapsed.count()*1./1000 << "s]\n";
    }

    TH1F tothist( "energy_total", "global MaGe energy spectrum", 2100, 0, 2.1 );
    TFile fileout("sumMaGe.root", "RECREATE");
    for ( auto& h : hist ) {
        h.Write();
        tothist.Add(&h);
    }
    tothist.Write();
    fileout.Close();

    return 0;
}
