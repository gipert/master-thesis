/* processData.cxx
 *
 * Saves energy spectra and live times on disk
 * using the class GERDA::DataReader
 *
 * Author: Luigi Pertoldi (luigi.pertoldi@pd.infn.it)
 * Created: 01/02/2017
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "TH1D.h"
#include "TFile.h"

#include "DataReader.h"

int main( int argc, char** argv ) {
    
    // copy arguments
    std::vector<std::string> args(argc);
    for ( int i = 0; i < argc; ++i ) args[i] = argv[i];

    // help
    if ( std::find(args.begin(), args.end(), "--help") != args.end() ) {
        std::ifstream helpFile( std::string(std::getenv("GERDACPTDIR")) + "/misc/help");
        if (helpFile.is_open()) std::cout << helpFile.rdbuf();
        else std::cerr << "help file not found!\n";
        return 0;
    }
/*
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
*/    
    std::string filename = std::string(std::getenv("GERDACPTDIR")) + "/data/sumData.root";

    std::string opt = "gauss";
    if ( std::find(args.begin(), args.end(), "--energy=zac") != args.end() or
         std::find(args.begin(), args.end(), "--energy=ZAC") != args.end() ) {
        opt = "zac";
        std::cout << "Using ZAC filter results...\n";
    }
 
 // ----------------------------------------------------------------------------------------------------
    // main reader object
    bool verbose = false;
    if ( std::find(args.begin(), args.end(), "--verbose") != args.end() ) verbose = true;
    GERDA::DataReader reader( std::string(std::getenv("GERDACPTDIR")) + "/misc/paths.txt", verbose, "GELATIO" );
    GERDA::DetectorSet set("GELATIO");

    // create output ROOT file and .txt file
    TFile file( filename.c_str(), "RECREATE" );
    filename.erase(filename.end()-4,filename.end());
    filename += "dat";
    std::ofstream textFile(filename.c_str());
    
    // retrieve energy spectrum
    std::vector<TH1D> energy;
    reader.CreateEnergyHist(opt);
    energy = reader.GetEnergyHist();

    TH1D energyBEGe("energyBEGeAll", "energyBEGeAll", 7500, 0, 7500);
    auto energyEnrCoax = reader.GetEnergyHistEnrCoax();

    for ( int i = 0; i < 40; ++i ) {
        // NOTE: excluding GTFs and GD02D
        if ( set.GetDetectorTypes()[i] == 1 and
             set.GetDetectorNames()[i] != "GD02D" ) {
            energyBEGe.Add(&energy[i]);
            energy[i].Write();
        }
        if ( set.GetDetectorTypes()[i] == 2 ) energy[i].Write();
    }
   
    // write on disk
    energyBEGe.Write();
    energyEnrCoax->Write();
    
    // retrieve time for each run [s]
    auto timeMap = reader.GetTimeMap();
    // write on disk
    for ( const auto& i : timeMap ) textFile << i.first << '\t' << i.second << '\n';
    
    textFile.close();

    std::string command = "telegram-send \"processData: completed\"";
    std::system(command.c_str());

    return 0;
}
