// DataReader.cxx
// 
// Author: Luigi Pertoldi
// Contact: luigi.pertoldi@pd.infn.it
// created: 31/01/2017

#include "DataReader.h"

#include <iostream>
#include <limits>
#include <memory>
#include <chrono>

// ROOT
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

// gerda-ADA
#include "FileMap.h"
#include "DataLoader.h"
#include "GETRunConfiguration.hh"

// other
#include "ProgressBar.h"

using namespace GERDA;

bool DataReader::kVerbosity = false;

DataReader::DataReader( std::string gerdaMetaPath, 
                        std::string gerdaDataPath,
                        std::string configListPath, 
                        std::string ordering ) :
    
    DetectorSet("GELATIO"),
    configList(configListPath.c_str()),
    kOrdering(ordering)

{ 
    if ( gerdaMetaPath.back() == '/' ) gerdaMetaPath.pop_back(); gerdaMetaDir = gerdaMetaPath;
    if ( gerdaDataPath.back() == '/' ) gerdaDataPath.pop_back(); gerdaDataDir = gerdaDataPath;

    energy.reserve(40);
    std::string histName;
    for ( int i = 0; i < 40; ++i ) {
        histName = "energy_";
        if ( detectorTypes[i] == 1 ) histName += "BEGe_";
        if ( detectorTypes[i] == 2 ) histName += "enrCoax_";
        if ( detectorTypes[i] == 3 ) histName += "natCoax_";
        histName += std::to_string(i);
        energy.emplace_back( histName.c_str(), histName.c_str(), 7500, 0, 7500 );
    }
    kMustResetEnergy = false;
}
// -------------------------------------------------------------------------------
DataReader::DataReader( std::string pathsFile, bool verbose, std::string ordering ) : 
    
    DetectorSet("GELATIO"),
    kOrdering(ordering)
{ 

    kVerbosity = verbose;

    std::ifstream input(pathsFile.c_str());
    if ( !input.is_open() ) std::cerr << "File with paths not found!\n";
    std::string metapath, datapath, configpath;
    input >> metapath >> datapath >> configpath;
    std::vector<unsigned int> runsToProcess;
    int value;
    while ( input >> value ) runsToProcess.push_back(value);
    input.close();

    configList.open(configpath.c_str());
    if ( metapath.back() == '/' ) metapath.pop_back(); gerdaMetaDir = metapath;
    if ( datapath.back() == '/' ) datapath.pop_back(); gerdaDataDir = datapath;

    energy.reserve(40);
    std::string histName;
    for ( int i = 0; i < 40; ++i ) {
        histName = "energy_";
        if ( detectorTypes[i] == 1 ) histName += "BEGe_";
        if ( detectorTypes[i] == 2 ) histName += "enrCoax_";
        if ( detectorTypes[i] == 3 ) histName += "natCoax_";
        histName += std::to_string(i);
        energy.emplace_back( histName.c_str(), histName.c_str(), 7500, 0, 7500 );
    }
    kMustResetEnergy = false;

    for ( auto& n : runsToProcess ) this->LoadRun(n);
    
        std::cout << "Loaded runs: ";
        for ( const auto& it : dataTreeMap ) std::cout << it.first << ' ';
        std::cout << ". Continue? [y/n, default=y] ";
        std::string line;
        std::getline( std::cin, line );
        if( !line.empty() and line != "y" ) { std::cout << "Aborting...\n"; throw -1; }
}
// -------------------------------------------------------------------------------
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
// -------------------------------------------------------------------------------
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

    if ( configFile.IsZombie() ) { std::cerr << "Run" << runID << ": config file not found!\n"; return false; }
    
    if (kVerbosity) std::cout << "Retrieving detector status...\n";
    std::unique_ptr<GETRunConfiguration> gtr(dynamic_cast<GETRunConfiguration*>(configFile.Get("RunConfiguration")));
       
    std::vector<int> detector_status( gtr->GetNDetectors(), 0 );
    for ( int i = 0; i < (int)detector_status.size(); ++i ) {
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

    /*auto tmp = loader.GetUniqueMasterChain();
    if ( tmp->IsZombie() or tmp == 0 ) {
        std::cerr << "Data chain of run" << runID << " is Zombie! Tree not loaded.\n";
        return false;
    }*/
    
    dataTreeMap.insert(std::make_pair( runID, std::unique_ptr<TChain>(loader.GetUniqueMasterChain()) )); 
    if (kVerbosity) std::cout << "Done.\n\n";

    return true;
}
// -------------------------------------------------------------------------------
void DataReader::CreateEnergyHist( std::string opt ) {

    if ( opt != "zac" and opt != "ZAC" and opt != "gauss" and opt != "GAUSS" ) {
        std::cout << opt << ": Unknown option, aborting...\n";
        return;
    }
   
    if (kMustResetEnergy) {
        std::cout << "Warning: the energy vector is non-empty, call DataReader::ResetEnergy. Aborting...\n"; 
        return;
    }
    
    int nTP;

    TTreeReader treereader;
    TTreeReaderValue<int> multiplicity     (treereader, "multiplicity.firedChannels");
    TTreeReaderValue<int> isTP             (treereader, "isTP.isTP");
    TTreeReaderValue<int> isVetoedInTime   (treereader, "isVetoedInTime.isvetoedintime");
    TTreeReaderArray<int> failedFlag       (treereader, "failedFlag");
    TTreeReaderArray<double> energyGauss   (treereader, "rawEnergyGauss");
    TTreeReaderArray<double> energyZAC     (treereader, "rawEnergyZAC");
    TTreeReaderArray<double> energyTot     (treereader, "energy");

    for ( const auto& it : dataTreeMap ) {
       
        if (kVerbosity) std::cout << "Initializing... " << std::flush;
        nTP = 0;
        auto& chain = it.second;
        treereader.SetTree(chain.get());

        ProgressBar bar(chain->GetEntries());
        std::cout << "processing run" << it.first << ": " << std::flush;
        
        auto start = std::chrono::system_clock::now();

        int i = 0;
        while (treereader.Next()) {
            
            bar.Update(i); i++;

            if (*isTP) nTP++;

            if ( !*isTP and !*isVetoedInTime and *multiplicity == 1 ) {
                for ( int det = 0; det < 40; det++ ) {

                    if ( opt == "gauss" or opt == "GAUSS" ) {
                    
                        // equivalent
                        /*if ( failedFlag->at(det) == 0 and 
                            detectorStatusMap[it.first][det] == 0 and
                            energyTot->at(det) > 0 and energyTot->at(det) < 10000 ) {
                            energy[det].Fill(energyGauss->at(det));
                        }*/
                    
                        if ( energyTot[det] > 0 and energyTot[det] < 10000 ) {
                            energy[det].Fill(energyTot[det]);
                        }
                    }

                    else if ( opt == "zac" or opt == "ZAC" ) {
                        
                        if ( failedFlag[det] == 0 and 
                            detectorStatusMap[it.first][det] == 0 and
                            energyTot[det] > 0 and energyTot[det] < 10000 ) {
                            energy[det].Fill(energyZAC[det]);
                        }

                    }
                }
            }
        }
        timeMap.insert(std::make_pair(it.first, nTP*20));

        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
        std::cout << " [" << elapsed.count()*1./1000 << "s]\n";
    }
    
    return;
}
// -------------------------------------------------------------------------------
void DataReader::ResetEnergy() {

    for ( auto& it : energy ) it.Reset();
    kMustResetEnergy = false;
    return;
}
// -------------------------------------------------------------------------------
std::vector<TH1D> DataReader::GetEnergyHist() {
    
    if ( kOrdering == "MaGeInput" ) {
        auto v = energy;
        GERDA::ReorderAsMaGe(v, "input");
        return v;
    }
    else if ( kOrdering == "MaGeOutput" ) {
        auto v = energy;
        GERDA::ReorderAsMaGe(v, "output");
        return v;
    }
    else return energy;
}
// -------------------------------------------------------------------------------
std::map<unsigned int, std::vector<int>> DataReader::GetDetectorStatusMap() { 
     
    if ( kOrdering == "MaGeInput" ) {
        auto v = detectorStatusMap;
        for ( auto& i : v ) GERDA::ReorderAsMaGe(i.second, "input");
        return v;
    }
    else if ( kOrdering == "MaGeOutput" ) {
        auto v = detectorStatusMap;
        for ( auto& i : v ) GERDA::ReorderAsMaGe(i.second, "output");
        return v;
    }
    else return detectorStatusMap;
}
// -------------------------------------------------------------------------------
unsigned int DataReader::GetTime() {
   
   unsigned int tmp = 0;
   for ( auto& it : timeMap ) tmp += it.second;
   return tmp;
}
// -------------------------------------------------------------------------------
std::unique_ptr<TH1D> DataReader::GetEnergyHistBEGe() const {
    
    std::unique_ptr<TH1D> tmp(new TH1D( "energyBEGeAll", "energyBegeAll", 7500, 0, 7500 ));
    if (energy.empty()) { std::cerr << "DataReader::CreateEnergyHist has not been called!\n"; return tmp; }
    
    for ( int i = 0; i < 40; ++i ) {
        if ( detectorTypes[i] == 1 ) tmp->Add(&energy[i]);
    }

    return tmp;
}
// -------------------------------------------------------------------------------
std::unique_ptr<TH1D> DataReader::GetEnergyHistEnrCoax() const {
    
    std::unique_ptr<TH1D> tmp(new TH1D( "energyEnrCoaxAll", "energyEnrCoaxAll", 7500, 0, 7500 ));
    if (energy.empty()) { std::cerr << "DataReader::CreateEnergyHist has not been called!\n"; return tmp; }
    
    for ( int i = 0; i < 40; ++i ) {
        if ( detectorTypes[i] == 2 ) tmp->Add(&energy[i]);
    }

    return tmp;
}
// -------------------------------------------------------------------------------
std::unique_ptr<TH1D> DataReader::GetEnergyHistNatCoax() const {
    
    std::unique_ptr<TH1D> tmp(new TH1D( "energyNatCoaxAll", "energyNatCoaxAll", 7500, 0, 7500 ));
    if (energy.empty()) { std::cerr << "DataReader::CreateEnergyHist has not been called!\n"; return tmp; }
    
    for ( int i = 0; i < 40; ++i ) {
        if ( detectorTypes[i] == 3 ) tmp->Add(&energy[i]);
    }

    return tmp;
}
// -------------------------------------------------------------------------------
TChain* DataReader::GetTreeFromRun( unsigned int runID ) const {

    auto result = dataTreeMap.find(runID);
    if ( result == dataTreeMap.end() ) {
        std::cout << "DataReader::GetTree: Run" << runID << " not loaded!\n";
        return nullptr;
    }

    return result->second.get();
}
// -------------------------------------------------------------------------------
TChain* DataReader::GetTree() {
    
    if (!dataTree) {
        if (kVerbosity) std::cout << "Creating Master Tree for the first time..." << std::endl;
        dataTree = std::unique_ptr<TChain>(new TChain());
        for ( auto& it : dataTreeMap ) dataTree->Add( it.second.get() );
    }

    return dataTree.get();
}
// -------------------------------------------------------------------------------
std::unique_ptr<TChain> DataReader::MoveTree() {

    if (!dataTree) this->GetTree();
    return std::move(dataTree);
}
// -------- end of class ------------------------------------------------------------
