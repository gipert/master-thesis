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

int main( int argc, char** argv ) {
    
    // copy arguments
    std::vector<std::string> args(argc);
    for ( int i = 0; i < argc; i++ ) args[i] = argv[i];

    // help
    if ( std::find(args.begin(), args.end(), "--help") != args.end() ) {
        std::ifstream helpFile("help");
        if (helpFile.is_open()) std::cout << helpFile.rdbuf();
        else std::cerr << "help file not found!\n";
        return 0;
    }
    
    // set verbose
    bool verbose = false;
    if ( std::find(args.begin(), args.end(), "--verbose") != args.end() ) verbose = true;

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

    // get configs
    std::ifstream input("paths.txt");
    if ( !input.is_open() ) { std::cerr << "File with paths not found! Aborting...\n"; return 0; }
    std::string metapath, datapath, configpath;
    input >> metapath >> datapath >> configpath;
    std::vector<unsigned int> runsToProcess;
    int value;
    while ( input >> value ) runsToProcess.push_back(value);
    
    // main reader object
    GERDA::DataReader reader( metapath, datapath, configpath );
    
    // load runs
    for ( auto& n : runsToProcess ) reader.LoadRun(n, verbose);
    
    // retrieve energy spectrum
    std::vector<TH1D> energy;
    reader.CreateEnergyHist();
    energy = reader.GetEnergyHist();

    TH1D* energyBEGe = reader.GetEnergyHistBEGe();
    TH1D* energyEnrCoax = reader.GetEnergyHistEnrCoax();
    TH1D* energyNatCoax = reader.GetEnergyHistNatCoax();

    for ( auto& i : runsToProcess ) std::cout << "Acquisition time for run" << i << ": " << reader.GetTimeHoursForRun(i) << " h\n";
    std::cout << "Total acquisition time: " << reader.GetTimeHours() << std::endl;
    std::cout << reader.GetTree()->GetEntries() << std::endl;
    reader.GetTree()->GetListOfFriends()->Print();

    std::string processedRunsStr;
    for ( auto& r : runsToProcess ) processedRunsStr += std::to_string(r) + ' ';
    TNamed processed_runs( "processed runs", processedRunsStr);

    TFile file( filename.c_str(), "RECREATE" );
    for ( const auto& it : energy ) it.Write();
    processed_runs.Write();
    energyBEGe->Write();
    energyEnrCoax->Write();
    energyNatCoax->Write();
    file.Close();
    
    delete energyBEGe;
    delete energyEnrCoax;
    delete energyNatCoax;

    return 0;
}
