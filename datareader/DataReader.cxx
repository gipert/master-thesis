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

// other
#include "progressbar.h"

using namespace GERDA;

bool DataReader::kVerbosity = false;

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
 
    energy.reserve(40);
    std::string histName;
    for ( int i = 0; i < 40; i++ ) {
        histName = "energy_";
        if ( detectorMatrix[i] == 1 ) histName += "BEGe_";
        if ( detectorMatrix[i] == 2 ) histName += "enrCoax_";
        if ( detectorMatrix[i] == 3 ) histName += "natCoax_";
        histName += std::to_string(i);
        energy.emplace_back( histName.c_str(), histName.c_str(), 7500, 0, 7500 );
    }
    kMustResetEnergy = false;
   
    dataTree = nullptr;
}

DataReader::~DataReader() {
    
    configList.close();
    for ( auto& ch : dataTreeMap ) delete ch.second;
    delete dataTree;
}

std::string DataReader::FindRunConfiguration( unsigned int runID ) {

    if ( !configList.is_open() ) return "filenotfound";
    configList.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    configList.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    
    std::string dummy;
    std::string result = "runnotregistered";
    while ( configList >> dummy ) {
        if ( (unsigned int)std::stoi(dummy) == runID ) {
            configList >> result;
            break;
        }
        configList.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    configList.clear();
    configList.seekg(0, std::ios::beg);
    return result;
}

bool DataReader::LoadRun( unsigned int runID ) {

    auto result = dataTreeMap.find(runID);
    if ( result != dataTreeMap.end() ) {
        std::cout << "DataReader::LoadRun: Run" << runID << " already loaded!\n";
        return false;
    }

    if (kVerbosity) {
        std::cout << "RUN" << runID << std::endl;
        std::cout << "Looking for config list file...\n";
    }
    std::string confName = FindRunConfiguration( runID );
    if ( confName == "filenotfound"     ) { std::cerr << "Config list file not found!\n"; return false; }
    if (kVerbosity) std::cout << "Looking for config file...\n";
    if ( confName == "runnotregistered" ) { std::cerr << "Run" << runID << ": runID not found in config list!\n"; return false; }

    std::string completePath = gerdaMetaDir + "/config/_aux/geruncfg/" + confName; 
    if (kVerbosity) std::cout << "Opening config file...\n";
    TFile configFile( completePath.c_str(), "READ" );
    
    //auto configFile = std::unique_ptr<TFile, decltype(&TFile::Close)>{ 
    //    TFile::Open(completePath.c_str,"READ"),
    //    &TFile::Close
    //};

    if ( configFile.IsZombie() ) { std::cerr << "Run" << runID << ": config file not found!\n"; return false; }
    
    if (kVerbosity) std::cout << "Retrieving detector status...\n";
    std::unique_ptr<GETRunConfiguration> gtr(dynamic_cast<GETRunConfiguration*>(configFile.Get("RunConfiguration")));
       
    std::vector<unsigned int> detector_status( gtr->GetNDetectors(), 0 );
    for ( int i = 0; i < (int)detector_status.size(); i++ ) {
        if      (  gtr->IsTrash(i) ) detector_status[i] = 2;
        else if ( !gtr->IsOn(i)    ) detector_status[i] = 1;
    }

    configFile.Close();
    detectorStatusMap.insert(std::make_pair( runID, detector_status ));
 
    if (kVerbosity) std::cout << "Looking for data files...\n";
    gada::FileMap myMap;
    myMap.SetRootDir(gerdaDataDir);
    std::string pathToListOfKeys = gerdaMetaDir + "/data-sets/phy/run00" + std::to_string(runID) + "-phy-analysis.txt";
    std::ifstream ftmp(pathToListOfKeys.c_str());
    if ( !ftmp.is_open() ) { std::cerr << pathToListOfKeys << " does not exist!\n"; return false; }
    myMap.BuildFromListOfKeys(pathToListOfKeys);

    if (kVerbosity) std::cout << "Loading trees...\n";
    gada::DataLoader loader;
    loader.AddFileMap(&myMap);
    if ( !loader.BuildTier3() ) {
        std::cerr << "DataLoader::BuildTier3 failed for run" << runID << ", tree not loaded.\n";
        return false;
    }

    if ( !loader.BuildTier4() ) {
        std::cerr << "DataLoader::BuildTier4 failed for run" << runID << ", tree not loaded.\n";
        return false;
    }

    auto tmp = loader.GetUniqueMasterChain();
    if ( tmp->IsZombie() or tmp == 0 ) {
        std::cerr << "Data chain of run" << runID << " is Zombie! Tree not loaded.\n";
        return false;
    }

    dataTreeMap.insert(std::make_pair( runID, tmp )); 
    if (kVerbosity) std::cout << "Done.\n\n";

    return true;
}

void DataReader::CreateEnergyHist() {
   
    if (kMustResetEnergy) {
        std::cout << "Warning: the energy vector is non-empty, call DataReader::ResetEnergy. Aborting...\n"; 
        return;
    }
    
    int nTP;
    int nEntries;
    int multiplicity, isTP, isVetoedInTime;
    std::vector<int>*    failedFlag = new std::vector<int>(40);
    std::vector<double>* energyGauss = new std::vector<double>(40);
    TChain* chain;

    for ( const auto& it : dataTreeMap ) {
        
        if (kVerbosity) std::cout << "Initialising... ";
        nTP = 0;
        chain = it.second;
        nEntries = chain->GetEntries();

        chain->SetBranchAddress("multiplicity"  , &multiplicity);
        chain->SetBranchAddress("rawEnergyGauss", &energyGauss);
        chain->SetBranchAddress("isTP"          , &isTP);
        chain->SetBranchAddress("isVetoedInTime", &isVetoedInTime);
        chain->SetBranchAddress("failedFlag"    , &failedFlag);

        ProgressBar bar(nEntries);
        std::cout << "processing run" << it.first << ": " << std::flush;
        bar.Init();
        
        for ( int e = 0; e < nEntries; e++ ) {
            
            bar.Update(e);
            chain->GetEntry(e);

            if (isTP) nTP++;

            if ( !isTP and !isVetoedInTime and multiplicity == 1 ) {
                for ( int det = 0; det < 40; det++ ) {
                    if ( !failedFlag->at(det) and detectorStatusMap[it.first][det] == 0 ) {
                        energy[det].Fill(energyGauss->at(det));
                    }
                }
            }
        }
        std::cout << std::endl;
        chain->ResetBranchAddresses();
        time.insert(std::make_pair(it.first, nTP*20));
    }
    
    delete failedFlag;
    delete energyGauss;

    return;
}

void DataReader::ResetEnergy() {

    for ( auto& it : energy ) it.Reset();
    kMustResetEnergy = false;
    return;
}

unsigned int DataReader::GetTime() {
   
   unsigned int tmp = 0;
   for ( auto& it : time ) tmp += it.second;
   return tmp;
}

TH1D* DataReader::GetEnergyHistBEGe() const {
    
    TH1D* tmp = new TH1D( "energyBEGeAll", "energyBegeAll", 7500, 0, 7500 );
    if (energy.empty()) { std::cerr << "DataReader::CreateEnergyHist has not been called!\n"; return tmp; }
    
    for ( int i = 0; i < 40; i++ ) {
        if ( detectorMatrix[i] == 1 ) tmp->Add(&energy[i]);
    }

    return tmp;
}

TH1D* DataReader::GetEnergyHistEnrCoax() const {
    
    TH1D* tmp = new TH1D( "energyEnrCoaxAll", "energyEnrCoaxAll", 7500, 0, 7500 );
    if (energy.empty()) { std::cerr << "DataReader::CreateEnergyHist has not been called!\n"; return tmp; }
    
    for ( int i = 0; i < 40; i++ ) {
        if ( detectorMatrix[i] == 2 ) tmp->Add(&energy[i]);
    }

    return tmp;
}

TH1D* DataReader::GetEnergyHistNatCoax() const {
    
    TH1D* tmp = new TH1D( "energyNatCoaxAll", "energyNatCoaxAll", 7500, 0, 7500 );
    if (energy.empty()) { std::cerr << "DataReader::CreateEnergyHist has not been called!\n"; return tmp; }
    
    for ( int i = 0; i < 40; i++ ) {
        if ( detectorMatrix[i] == 3 ) tmp->Add(&energy[i]);
    }

    return tmp;
}


TChain* DataReader::GetTreeFromRun( unsigned int runID ) const {

    auto result = dataTreeMap.find(runID);
    if ( result == dataTreeMap.end() ) {
        std::cout << "DataReader::GetTree: Run" << runID << " not loaded!\n";
        return nullptr;
    }

    return result->second;
}

TChain* DataReader::GetTree() {
    
    if (!dataTree) {
        if (kVerbosity) std::cout << "Creating Master Tree for the first time..." << std::endl;
        TChain* chain = new TChain();
        for ( auto& it : dataTreeMap ) chain->Add( it.second );
        dataTree = chain;
    }

    return dataTree;
}

/*TChain* DataReader::GetUniqueTree() const {
    TChain chain( *GetTree() );
    return &chain;
}*/
