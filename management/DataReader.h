// DataReader.h
// 
// class to load and read GERDA datasets 
// written on top of gerda-ADA classes
// 
// Author: Luigi Pertoldi
// Contact: luigi.pertoldi@pd.infn.it
// created: 31/01/2017
//
// NOTES: - the internal DetectorSet has always "GELATIO" odering,
//          also if DataReader is "MaGe", keep in mind. Decided to 
//          make inheritance private.

#ifndef DATA_READER__
#define DATA_READER__

#include <map>
#include <string>
#include <fstream>
#include <memory>

// ROOT
#include "TChain.h"
#include "TH1D.h"

// other
#include "DetectorSet.h"

namespace GERDA {

  class DataReader : private DetectorSet {
    
      public:
      
      // delete default constructor
      DataReader()                             = delete;
      // delete copy constructor/assignement
      DataReader           (DataReader const&) = delete;
      DataReader& operator=(DataReader const&) = delete;
      // default move constructor/assignement, redundand
      DataReader           (DataReader&&)      = default;
      DataReader& operator=(DataReader&&)      = default;
     
      // set paths, need then to call LoadRun to complete the configuration
      DataReader( std::string gerdaMetaPath,    // location of gerda-metadata repo
                  std::string gerdaDataPath,    // location of gerda-data folder
                  std::string configListPath,   // location of runconfiguration_mod.db
                  std::string ordering = "GELATIO" );       // "GELATIO" or "MaGeInput" or "MaGeOutput"
      // shortcut: set paths and load runs stored in an input file
      DataReader( std::string pathsFile, bool verbose = false, std::string ordering = "GELATIO" );
      
      // override default destructor
      ~DataReader() { configList.close(); }

      // load tree in DataTreeMap (optional: verbose mode)
      bool LoadRun( unsigned int runID );
      // get energy histogram from all runs with default cuts:
      int CreateEnergyHist( std::string opt = "gauss" );
      // clear energy histograms
      void ResetEnergy();
      // get vector with energy spectra for each detector
      std::vector<TH1D> GetEnergyHist();
      // get detector status
      std::map<unsigned int, std::vector<int>> GetDetectorStatusMap();
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
      // get non-owning pointers to trees
      // WARNING: deleted when the DataReader object goes out of scope
      TChain* GetTreeFromRun( unsigned int runID ) const;
      TChain* GetTree();
      // get owning pointer to global tree
      std::unique_ptr<TChain> MoveTree();
      
      static bool kVerbosity;

      private:
    
      // detector status: 0 = ON, 1 = AC, 2 = OFF
      std::map<unsigned int, std::vector<int>> detectorStatusMap;
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
      // vector with energy histograms for each detector
      // filled by GetEnergyHist();
      std::vector<TH1D> energy;
      // time in [s]
      std::map<unsigned int, unsigned int> timeMap;
 
      // flag for multiple calling of CreateEnergyhist
      bool kMustResetEnergy;
      const std::string kOrdering;
   
      // find run configuration in the config list file
      std::string FindRunConfiguration( unsigned int runID );
  }; 
}

#endif
