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
#include "TNamed.h"
#include "TParameter.h"

int main( int argc, char** argv ) {
    
    // copy arguments
    std::vector<std::string> args(argc);
    for ( int i = 0; i < argc; i++ ) args[i] = argv[i];

    // help
    if ( std::find(args.begin(), args.end(), "--help") != args.end() ) {
        std::ifstream helpFile("misc/help");
        if (helpFile.is_open()) std::cout << helpFile.rdbuf();
        else std::cerr << "help file not found!\n";
        return 0;
    }

    // get output filename
    std::string filename;
    for ( auto& str : args ) {
        if ( str.find(".root") != std::string::npos ) { 
            filename = str;
            break;
        }
    }

    if ( filename.empty() ) {
        std::cerr << "Please provide a valid .root output filename.\n";
        return -1;
    }

    std::string opt = "gauss";
    if ( std::find(args.begin(), args.end(), "--energy=zac") != args.end() or
         std::find(args.begin(), args.end(), "--energy=ZAC") != args.end() ) {
        opt = "zac";
        std::cout << "Using ZAC filter results...\n";
    }
    
    // get configs
    std::ifstream input("misc/paths.txt");
    if ( !input.is_open() ) { std::cerr << "File with paths not found! Aborting...\n"; return 0; }
    std::string metapath, datapath, configpath;
    input >> metapath >> datapath >> configpath;
    std::vector<unsigned int> runsToProcess;
    int value;
    while ( input >> value ) runsToProcess.push_back(value);
    input.close();
    
    // main reader object
    GERDA::DataReader reader( metapath, datapath, configpath );    
    // set verbose
    if ( std::find(args.begin(), args.end(), "--verbose") != args.end() ) GERDA::DataReader::kVerbosity = true;
    
    // load runs
    for ( auto& n : runsToProcess ) reader.LoadRun(n);
    
    // retrieve energy spectrum
    std::vector<TH1D> energy;
    reader.CreateEnergyHist(opt);
    energy = reader.GetEnergyHist();

    auto energyBEGe = reader.GetEnergyHistBEGe();
    auto energyEnrCoax = reader.GetEnergyHistEnrCoax();
    auto energyNatCoax = reader.GetEnergyHistNatCoax();
   
    TParameter<float> time( "total_acq_time_in_h", reader.GetTimeHours() );

    std::string processedRunsStr;
    for ( auto& r : runsToProcess ) processedRunsStr += std::to_string(r) + ' ';
    TNamed processed_runs( "processed runs", processedRunsStr);
    
    TFile file( filename.c_str(), "RECREATE" );
    for ( const auto& it : energy ) it.Write();
    processed_runs.Write();
    time.Write();
    energyBEGe->Write();
    energyEnrCoax->Write();
    energyNatCoax->Write();
    file.Close();
    
    return 0;
}

    /* NOT WORKING
    // retrieve run IDs from command line
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


