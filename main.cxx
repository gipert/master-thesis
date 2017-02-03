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

    reader.LoadRun(53, verbose);

    auto chain = reader.GetGlobalTree();
    chain->Print("all");

    
    
    return 0;
}
