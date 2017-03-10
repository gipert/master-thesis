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

      DataReader( std::string gerdaMetaPath,    // location of gerda-metadata repo
                  std::string gerdaDataPath,    // location of gerda-data folder
                  std::string configListPath ); // location of runconfiguration_mod.db
      
      // override default destructor
      ~DataReader();

      // load tree in DataTreeMap (optional: verbose mode)
      bool LoadRun( unsigned int runID );
      // get energy histogram from all runs with default cuts:
      void CreateEnergyHist( std::string opt = "gauss" );
      // clear energy histograms
      void ResetEnergy();
      // get non-owned vector with energy spectra for each detector
      std::vector<TH1D> GetEnergyHist() const { return energy; }
      // get acquisition time (in minutes)
      unsigned int GetTimeForRun( unsigned int runID ) { return time.at(runID); }
      unsigned int GetTime();
      float GetTimeHoursForRun( unsigned int runID ) { return (this->GetTimeForRun(runID))*1./3600; }
      float GetTimeHours() { return (this->GetTime())*1./3600; }
      // get owning pointers for histograms: 
      // WARNING: delete them to prevent memory leaks
      TH1D* GetEnergyHistBEGe()    const;
      TH1D* GetEnergyHistEnrCoax() const;
      TH1D* GetEnergyHistNatCoax() const;
      TH1D* GetEnergyHistAll()     const;
      // get non-owning pointers to trees
      // WARNING: deleted when the DataReader object goes out of scope
      TChain* GetTreeFromRun( unsigned int runID ) const;
      TChain* GetTree();
      // TODO: TChain* GetUniqueTree() const;
      
      static bool kVerbosity;

      private:

      // description of detectors types in detector strings
      std::vector<unsigned int> detectorMatrix;
      // config list file
      std::ifstream configList;
      // paths to gerda-metadata repo and data directory
      // note: the data directory structure must match the 
      // one requested by gerda-ADA
      std::string gerdaMetaDir;
      std::string gerdaDataDir;
      // map with trees, the key is the run ID
      std::map<unsigned int, TChain*> dataTreeMap;
      TChain* dataTree;
      // map with the detector's status, the key is the run ID
      std::map<unsigned int, std::vector<unsigned int>> detectorStatusMap;
      // vector with energy histograms for each detector
      // filled by GetEnergyHist();
      std::vector<TH1D> energy;
      // flag for multiple calling of CreateEnergyhist
      bool kMustResetEnergy;
      std::map<unsigned int, unsigned int> time;
    
      // find run configuration in the config list file
      std::string FindRunConfiguration( unsigned int runID );
  };
}

#endif
