/* DetectorSet.h
 *
 * class to store relevant informations concearning a set of
 * detectors and methods to calculate useful quantities
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 22/03/2017
 *
 * Ordering convention: GELATIO scheme
 *
 */

#ifndef DETECTOR_SET__
#define DETECTOR_SET__

#include <vector>
#include <string>

namespace GERDA {

  class DetectorSet {

      public:

      DetectorSet( std::string naming = "GELATIO" );
      ~DetectorSet()                             = default;
      // delete copy constructor/assignement
      DetectorSet           (DetectorSet const&) = delete;
      DetectorSet& operator=(DetectorSet const&) = delete;
      // default move constructor/assignement, redundand
      DetectorSet           (DetectorSet&&)      = default;
      DetectorSet& operator=(DetectorSet&&)      = default;

      
      std::vector<int>   GetDetectorTypes() const { return detectorTypes; }
      // get mass [g]
      std::vector<int>   GetMass()          const { return mass; }
      // get total volume, active volume, dead volume [cm^3]
      std::vector<float> GetVolume()        const;
      std::vector<float> GetActiveVolume()  const;
      std::vector<float> GetDeadVolume()    const;
      // get number of 76Ge isotopes:
      std::vector<double> GetN76Ge()        const;

      protected:

      // detectors types in detector strings:
      // 1 = BEGe, 2 = enrCOAX, 3 = natCOAX
      std::vector<int> detectorTypes;
      // masses
      std::vector<int> mass; // [g]
      // active volume fractions
      std::vector<float> fractionAV;
      // 76Ge enrichment fractions
      std::vector<float> fractionEnr;
      // mean molar mass
      std::vector<float> meanMolarMass;
      // mass densities
      const float natGeDensity = 5.32; // [g/cm^3]
      const float enrGeDensity = 5.54;
      // Avogadro number
      const double NAv = 6.02214E23;
  };
  
  // order vectors in the MaGe naming scheme
  template<typename T>
  void ReorderAsMaGe( std::vector<T>& v ) {
    
    int c = 0;
    auto v_ = v;
    for ( int i = 37; i <= 39 ; i++ ) { v[c] = v_[i]; c++; }
    for ( int i = 8 ; i <= 10 ; i++ ) { v[c] = v_[i]; c++; }
    for ( int i = 27; i <= 29 ; i++ ) { v[c] = v_[i]; c++; }
    v[c] = v_[36]; c++;
    for ( int i = 0 ; i <= 7  ; i++ ) { v[c] = v_[i]; c++; }
    for ( int i = 11; i <= 26 ; i++ ) { v[c] = v_[i]; c++; }
    for ( int i = 30; i <= 35 ; i++ ) { v[c] = v_[i]; c++; }
  }
}

#endif
