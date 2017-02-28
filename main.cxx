// test main for GERDA::DataReader class
//
// Author: Luigi Pertoldi (luigi.pertoldi@pd.infn.it)
// Created: 01/02/2017

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "DataReader.h"
#include "TChain.h"
#include "TH1D.h"
#include "TFile.h"

int main( int argc, char** argv ) {
    
    // copy arguments
    std::vector<std::string> args(argc);
    for ( int i = 0; i < argc; i++ ) args[i] = argv[i];
    
    // set verbose
    bool verbose = false;
    if ( std::find( args.begin(), args.end(), "--verbose") != args.end() ) verbose = true;
    
    /* NOT WORKING
    // retrieve run IDs
    std::vector<int> runsToProcess;
    std::string tmp;
    for ( auto& it : args ) {
        if ( it.front() == '{' and it.back() == '}' ) {
            for ( int i = 1; i < it.size()-1; i++ ) {
                if ( it[i] != ',' ) tmp.append(it[i]);
                if ( it[i] == ',' ) {
                    if ( tmp.find('-') != std::string::npos ) {
                        std::string tmp1, tmp2;
                        tmp1.append(tmp[0]);
                        tmp1.append(tmp[1]);
                        tmp2.append(tmp[3]);
                        tmp2.append(tmp[4]);
                        for ( int j = std::stoi(tmp1); j <= std::stoi(tmp2); j++ ) {
                            runsToProcess.push_back(j);
                        }
                        tmp1.clear();
                        tmp2.clear();
                    }
                    else {
                        runsToProcess.push_back(std::stoi(tmp));
                        tmp.clear();
                    }
                }
            }
        }

        else { std::cout << "bad input! See --help"<< std::endl; return -1; }
    }
    for ( auto& i : runsToProcess ) std::cout << i << " ";
    */

    std::ifstream input("paths.txt");
    if ( !input.is_open() ) { std::cout << "File with paths not found! Aborting...\n"; return 0; }
    
    std::string metapath, datapath, configpath;
    input >> metapath >> datapath >> configpath;

    std::vector<int> runsToProcess;
    int value;
    while ( input >> value ) runsToProcess.push_back(value);
    
    GERDA::DataReader reader( metapath, datapath, configpath);
    
    for ( auto& it : runsToProcess ) reader.LoadRun(it, verbose);

    std::vector<TH1D> energy;
    energy = reader.GetEnergyHist();

    TFile file( "results.root", "RECREATE" );
    for ( const auto& it : energy ) it.Write();
    file.Close();
    
    return 0;
}
