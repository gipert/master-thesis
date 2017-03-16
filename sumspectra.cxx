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
    auto maxvolumeAV = std::max_element(AV.begin(), AV.end());
    auto maxvolumeDV = std::max_element(DV.begin(), DV.end());
    float maxvolume = *maxvolumeAV > *maxvolumeDV ? *maxvolumeAV : *maxvolumeDV;

    // define reading objects
    std::unique_ptr<TFile> file;
    TTreeReader treereader("fTree", file.get());
    
    TTreeReaderArray<int>   det_id(treereader, "det_id");
    TTreeReaderArray<float> det_edep(treereader, "det_edep");

    // final histogram
    std::vector<TH1F> hist;
    for ( int i = 0; i < 40; i++ ) {
        hist.emplace_back(Form("energy_det_id%i",i), "global MaGe spectra", 7500, 0, 7500);
    }
    
    ProgressBar bar;
    std::string filename;
    int nentries;
    int size;
    
    // loop over enrCOAX files
    for ( int i = 4; i <= 10; i++ ) {

        // active volume
        filename = "2nbbLV_AV_det11_" + std::to_string(i) + ".root";
        file = std::unique_ptr<TFile>{ TFile::Open(filename.c_str(), "READ") };
        
        // fill with correct entries
        nentries = treereader.GetEntries(true)*(1 - AV[i-1]/maxvolume);
        bar.SetNIter(nentries);
        bar.Init();
        for ( int j = 0; j < nentries; j++ ) {
            bar.Update(j);
            treereader.Next();
            size = det_id.GetSize();
            for ( int k = 0; k < size; k++ ) hist[det_id[k]].Fill(det_edep[k]);
        }

        // dead volume
        filename = "2nbbLV_DV_det11_" + std::to_string(i) + ".root";
        file = std::unique_ptr<TFile>{ TFile::Open(filename.c_str(), "READ") };
        
        // fill with correct entries
        nentries = treereader.GetEntries(true)*(1 - DV[i-1]/maxvolume);
        bar.SetNIter(nentries);
        bar.Init();
        for ( int j = 0; j < nentries; j++ ) {
            bar.Update(j);
            treereader.Next();
            size = det_id.GetSize();
            for ( int k = 0; k < size; k++ ) hist[det_id[k]].Fill(det_edep[k]);
        }
    }

    // loop over BEGe
    for ( int i = 1; i <= 30; i++ ) {

        // active volume
        filename = "2nbbLV_AV_det5_" + std::to_string(i) + ".root";
        file = std::unique_ptr<TFile>{ TFile::Open(filename.c_str(), "READ") };
        
        // fill with correct entries
        nentries = treereader.GetEntries(true)*(1 - AV[i+9]/maxvolume);
        bar.SetNIter(nentries);
        bar.Init();
        for ( int j = 0; j < nentries; j++ ) {
            bar.Update(j);
            treereader.Next();
            size = det_id.GetSize();
            for ( int k = 0; k < size; k++ ) hist[det_id[k]].Fill(det_edep[k]);
        }

        // dead volume
        filename = "2nbbLV_DV_det5_" + std::to_string(i) + ".root";
        file = std::unique_ptr<TFile>{ TFile::Open(filename.c_str(), "READ") };
        
        // fill with correct entries
        nentries = treereader.GetEntries(true)*(1 - DV[i+9]/maxvolume);
        bar.SetNIter(nentries);
        bar.Init();
        for ( int j = 0; j < nentries; j++ ) {
            bar.Update(j);
            treereader.Next();
            size = det_id.GetSize();
            for ( int k = 0; k < size; k++ ) hist[det_id[k]].Fill(det_edep[k]);
        }
    }
    
    TFile fileout("sumMaGe.root", "RECREATE");
    for ( auto& h : hist ) h.Write();
    fileout.Close();

    return 0;
}
