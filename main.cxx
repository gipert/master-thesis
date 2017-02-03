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

    for ( int i = 53; i < 66; i++) reader.LoadRun(i, verbose);
    
    std::cout << "Creating global tree...\n" << std::flush;
    auto chain = reader.GetGlobalTree();
    std::cout << "Getting number of entries...\n" << std::flush;
    std::cout << chain->GetEntries() << std::endl;

    
    
    return 0;
}
