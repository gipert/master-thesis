// DataReader.cxx
// 
// Author: Luigi Pertoldi
// Contact: luigi.pertoldi@pd.infn.it
// created: 31/01/2017

#include "DataReader.h"

#include <iostream>
#include <limits>
#include <memory>

// ROOT
#include "TFile.h"

// gerda-ADA
#include "FileMap.h"
#include "DataLoader.h"
#include "GETRunConfiguration.hh"

using namespace GERDA;

DataReader::DataReader( std::string gerdaMetaPath, 
                        std::string gerdaDataPath,
                        std::string configListPath ) :
    
    configList(configListPath.c_str())
{ 
    if ( gerdaMetaPath.back() == '/' ) gerdaMetaPath.pop_back(); gerdaMetaDir = gerdaMetaPath;
    if ( gerdaDataPath.back() == '/' ) gerdaDataPath.pop_back(); gerdaDataDir = gerdaDataPath;
    
    detectorMatrix = { 1,1,1,1,1,1,1,1, /*string1*/
                       2,2,2,           /*string2*/
                       1,1,1,1,1,1,1,1, /*string3*/
                       1,1,1,1,1,1,1,1, /*string4*/
                       2,2,2,           /*string5*/
                       1,1,1,1,1,1,2,   /*string6*/
                       3,3,3            /*string7*/ }; // 1 BEGe, 2 enrCoax, 3 natCoax
}

DataReader::~DataReader() {
    
    configList.close();
    for ( auto& ch : dataTree ) delete ch.second;  
}

std::string DataReader::FindRunConfiguration( int runID ) {

    if ( !configList.is_open() ) return "filenotfound";
    configList.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    configList.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    
    std::string dummy;
    std::string result = "runnotregistered";
    while ( configList >> dummy ) {
        if ( std::stoi(dummy) == runID ) {
            configList >> result;
            break;
        }
        configList.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    configList.clear();
    configList.seekg(0, std::ios::beg);
    return result;
}

bool DataReader::LoadRun( int runID , bool verbose ) {

    if ( verbose == true ) {
        std::cout << "RUN" << runID << std::endl;
        std::cout << "Looking for config list file...\n";
    }
    std::string confName = FindRunConfiguration( runID );
    if ( confName == "filenotfound"     ) { std::cout << "Config list file not found!\n"; return false; }
    if ( verbose == true ) std::cout << "Looking for config file...\n";
    if ( confName == "runnotregistered" ) { std::cout << "Run" << runID << ": runID not found in config list!\n"; return false; }

    std::string completePath = gerdaMetaDir + "/config/_aux/geruncfg/" + confName; 
    if ( verbose == true ) std::cout << "Opening config file...\n";
    TFile configFile( completePath.c_str(), "READ" );
    
    //auto configFile = std::unique_ptr<TFile, decltype(&TFile::Close)>{ 
    //    TFile::Open(completePath.c_str,"READ"),
    //    &TFile::Close
    //};

    if ( configFile.IsZombie() ) { std::cout << "Run" << runID << ": config file not found!\n"; return false; }
    
    if ( verbose == true ) std::cout << "Retrieving detector status...\n";
    std::unique_ptr<GETRunConfiguration> gtr(dynamic_cast<GETRunConfiguration*>(configFile.Get("RunConfiguration")));
       
    std::vector<int> detector_status( gtr->GetNDetectors(), 0 );
    for ( int i = 0; i < (int)detector_status.size(); i++ ) {
        if ( gtr->IsTrash(i) ) detector_status[i] = 2;
        if ( gtr->IsOn(i)    ) detector_status[i] = 1;
    }
    configFile.Close();
 
    if ( verbose == true ) std::cout << "Looking for data files...\n";
    gada::FileMap myMap;
    myMap.SetRootDir(gerdaDataDir);
    std::string pathToListOfKeys = gerdaMetaDir + "/data-sets/phy/run00" + std::to_string(runID) + "-phy-analysis.txt";
    std::ifstream ftmp(pathToListOfKeys.c_str());
    if ( !ftmp.is_open() ) { std::cout << pathToListOfKeys << " does not exist!\n"; return false; }
    myMap.BuildFromListOfKeys(pathToListOfKeys);

    if ( verbose == true ) std::cout << "Loading trees...\n";
    gada::DataLoader loader;
    loader.AddFileMap(&myMap);
    if ( !loader.BuildTier3() ) {
        std::cout << "DataLoader::BuildTier3 failed for run" << runID << ", tree not loaded.\n";
        return false;
    }

    if ( !loader.BuildTier4() ) {
        std::cout << "DataLoader::BuildTier4 failed for run" << runID << ", tree not loaded.\n";
        return false;
    }

    dataTree.insert(std::make_pair( runID, loader.GetUniqueMasterChain() )); 
    if ( verbose == true ) std::cout << "Done.\n\n";

    return true;
}

TChain* DataReader::GetTree( int runID ) const {

    auto result = dataTree.find(runID);
    if ( result == dataTree.end() ) {
        std::cout << "DataReader::GetTree: Run" << runID << " not loaded!\n";
        return nullptr;
    }

    return result->second;
}

TChain* DataReader::GetGlobalTree() const {

    TChain* chain = new TChain;
    for ( auto& it : dataTree ) chain->AddFriend( it.second );
    return chain;
}
