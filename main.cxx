// test main for GERDA::DataReader class
//
// Author: Luigi Pertoldi (luigi.pertoldi@pd.infn.it)
// Created: 01\02/2017

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "DataReader.h"
#include "TChain.h"
#include "TH1D.h"
#include "TFile.h"

int main( int argc, char** argv ) {

    std::vector<std::string> args(argc);
    for ( int i = 0; i < argc; i++ ) args[i] = argv[i];

    bool verbose = false;
    if ( argc > 1 and args[1] == "--verbose" ) verbose = true;

    std::ifstream input("paths.txt");
    if ( !input.is_open() ) { std::cout << "File with paths not found! Aborting...\n"; return 0; }
    
    std::string metapath, datapath, configpath;
    input >> metapath >> datapath >> configpath;
    
    GERDA::DataReader reader( metapath, datapath, configpath);

    reader.LoadRun(53, verbose);
    
    std::cout << "Getting tree...\n" << std::flush;
    auto chain = reader.GetTreeFromRun(53);
    std::cout << "Getting number of entries...\n" << std::flush;
    std::cout << chain->GetEntries() << std::endl;

    std::vector<TH1D> energy;
    energy = reader.GetEnergyHist();

    TFile file( "results.root", "RECREATE" );
    for ( const auto& it : energy ) it.Write();
    file.Close();
    
    return 0;
}
