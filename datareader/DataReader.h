// DataReader.h
// 
// class to load and read GERDA datasets 
// written on top of gerda-ADA classes
// 
// Author: Luigi Pertoldi
// Contact: luigi.pertoldi@pd.infn.it
// created: 31/01/2017

#ifndef DATA_READER__
#define DATA_READER__

#include <map>
#include <string>
#include <fstream>
#include <memory>

// ROOT
#include "TChain.h"
#include "TH1D.h"

namespace GERDA {

  class DataReader {
    
      // delete default constructor
      DataReader()                             = delete;
      // delete copy constructor/assignement
      DataReader           (DataReader const&) = delete;
      DataReader& operator=(DataReader const&) = delete;
      // default move constructor/assignement, redundand
      DataReader           (DataReader&&)      = default;
      DataReader& operator=(DataReader&&)      = default;

      public:
      
      // set paths, need then to call LoadRun to complete the configuration
      DataReader( std::string gerdaMetaPath,    // location of gerda-metadata repo
                  std::string gerdaDataPath,    // location of gerda-data folder
                  std::string configListPath ); // location of runconfiguration_mod.db
      // shortcut: set paths and load runs stored in an input file
      DataReader( std::string pathsFile, bool verbose = false );
      
      // override default destructor
      ~DataReader() { configList.close(); }

      // load tree in DataTreeMap (optional: verbose mode)
      bool LoadRun( unsigned int runID );
      // get energy histogram from all runs with default cuts:
      void CreateEnergyHist( std::string opt = "gauss" );
      // clear energy histograms
      void ResetEnergy();
      // get vector with energy spectra for each detector
      std::vector<TH1D> GetEnergyHist() const { return energy; }
      // get detector status map
      std::map<unsigned int, std::vector<unsigned int>> GetDetectorStatusMap() const { return detectorStatusMap; }
      // get acquisition time (in minutes)
      unsigned int GetTimeForRun( unsigned int runID ) { return timeMap.at(runID); }
      unsigned int GetTime();
      std::map<unsigned int, unsigned int> GetTimeMap() const { return timeMap; };
      float GetTimeHoursForRun( unsigned int runID ) { return (this->GetTimeForRun(runID))*1./3600; }
      float GetTimeHours() { return (this->GetTime())*1./3600; }
      // get owning pointers for energy histograms:
      std::unique_ptr<TH1D> GetEnergyHistBEGe()    const;
      std::unique_ptr<TH1D> GetEnergyHistEnrCoax() const;
      std::unique_ptr<TH1D> GetEnergyHistNatCoax() const;
      std::unique_ptr<TH1D> GetEnergyHistAll()     const;
      // get total volume, active volume, dead volume cm^3
      std::vector<float> GetVolume( std::string opt = "" ) const;
      std::vector<float> GetActiveVolume( std::string opt = "" ) const;
      std::vector<float> GetDeadVolume( std::string opt = "" ) const;
      std::vector<int>   GetMass( std::string opt = "" ) const;
      // get non-owning pointers to trees
      // WARNING: deleted when the DataReader object goes out of scope
      TChain* GetTreeFromRun( unsigned int runID ) const;
      TChain* GetTree();
      // get owning pointer to global tree
      std::unique_ptr<TChain> MoveTree();
      
      static bool kVerbosity;

      private:

      // description of detectors types in detector strings
      const std::vector<unsigned int> detectorMatrix;
      // total mass and active volume fraction:
      const std::vector<int>   mass; // g
      const std::vector<float> fractionAV;
      const float natGeDensity = 5.32; // g/cm^3
      const float enrGeDensity = 5.54;
      // config list file
      std::ifstream configList;
      // paths to gerda-metadata repo and data directory
      // note: the data directory structure must match the 
      // one requested by gerda-ADA
      std::string gerdaMetaDir;
      std::string gerdaDataDir;
      // map with trees, the key is the run ID
      std::map<unsigned int, std::unique_ptr<TChain>> dataTreeMap;
      std::unique_ptr<TChain> dataTree;
      // map with the detector's status, the key is the run ID
      // 0 = ON, 1 = AC, 2 = OFF
      std::map<unsigned int, std::vector<unsigned int>> detectorStatusMap;
      // vector with energy histograms for each detector
      // filled by GetEnergyHist();
      std::vector<TH1D> energy;
      // flag for multiple calling of CreateEnergyhist
      bool kMustResetEnergy;
      std::map<unsigned int, unsigned int> timeMap;
    
      // find run configuration in the config list file
      std::string FindRunConfiguration( unsigned int runID );
  };
  
  // order vectors in the MaGe input naming scheme
  template<typename T>
  void ReorderAsMaGeInput( std::vector<T>& v ) {
    
    // reorder as MaGe inpunt naming convention
    int c = 0;
    auto v_ = v;
    for ( int i = 37; i <= 39 ; i++ ) { v[c] = v_[i]; c++; }
    for ( int i = 8 ; i <= 10 ; i++ ) { v[c] = v_[i]; c++; }
    for ( int i = 27; i <= 29 ; i++ ) { v[c] = v_[i]; c++; }
    v[c] = v_[36]; c++;
    for ( int i = 0 ; i <= 7  ; i++ ) { v[c] = v_[i]; c++; }
    for ( int i = 11; i <= 26 ; i++ ) { v[c] = v_[i]; c++; }
    for ( int i = 30; i <= 35 ; i++ ) { v[c] = v_[i]; c++; }

    return;
  }
}

#endif
