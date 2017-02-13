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

      // load tree in dataTree (optional: verbose mode)
      bool LoadRun( int runID , bool verbose = false );
      // get energy histogram with default cuts:
      std::vector<TH1D> GetEnergyHist();
      // retrieve tree
      // WARNING: deleted when the DataReader object goes out of scope
      TChain* GetTreeFromRun( int runID ) const;
      TChain* GetTree() const;
      // TODO: GetOwningTree

      private:

      // description of detectors types in detector strings
      std::vector<int> detectorMatrix;
      // config list file
      std::ifstream configList;
      // paths to gerda-metadata repo and data directory
      // note: the data directory structure must match the 
      // one requested by gerda-ADA
      std::string gerdaMetaDir;
      std::string gerdaDataDir;
      // map with trees, the key is the run ID
      std::map<int, TChain*> dataTree;
      // map with detector's status, the key is the run ID
      std::map<int, std::vector<int>> detectorStatusMap;
    
      // find run configuration in the config list file
      std::string FindRunConfiguration( int runID );
  };
}

#endif
