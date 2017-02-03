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
    
    return result;
}

bool DataReader::LoadRun( int runID ) {

    std::string confName = FindRunConfiguration( runID );
    if ( confName == "filenotfound"     ) { std::cout << "Config list file not found!\n"; return false; }
    if ( confName == "runnotregistered" ) { std::cout << "Run" << runID << ": runID not found in config list!\n"; return false; }

    std::string completePath = gerdaMetaDir + "/config/_aux/geruncfg/" + confName;
    //auto configFile = std::unique_ptr<TFile, decltype(&TFile::Close)>{ 
    //    TFile::Open(completePath.c_str,"READ"),
    //    &TFile::Close
    //};

    TFile configFile( completePath.c_str(), "READ" );
    if ( configFile.IsZombie() ) { std::cout << "Run" << runID << ": config file not found!\n"; return false; }
    
    std::unique_ptr<GETRunConfiguration> gtr(dynamic_cast<GETRunConfiguration*>(configFile.Clone("RunConfiguration")));

    std::vector<int> detector_status( gtr->GetNDetectors(), 0 );
    for ( auto i : detector_status ) {
        if ( gtr->IsTrash(i) ) detector_status[i] = 2;
        if ( gtr->IsOn(i)    ) detector_status[i] = 1;
    }

    configFile.Close();
    return true;
}
