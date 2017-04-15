/* DetectorSet.h
 *
 * class to store relevant informations concearning a set of
 * detectors and methods to calculate useful quantities
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 22/03/2017
 *
 */

#ifndef DETECTOR_SET__
#define DETECTOR_SET__

#include <iostream>
#include <vector>
#include <string>

namespace GERDA {

  class DetectorSet {

      public:

      DetectorSet( std::string naming );
      DetectorSet()                              = delete;
      ~DetectorSet()                             = default;
      // delete copy constructor/assignement
      DetectorSet           (DetectorSet const&) = delete;
      DetectorSet& operator=(DetectorSet const&) = delete;
      // default move constructor/assignement, redundand
      DetectorSet           (DetectorSet&&)      = default;
      DetectorSet& operator=(DetectorSet&&)      = default;

      
      std::vector<std::string> GetDetectorNames()   const { return names; }
      // 1 = BEGe, 2 = enrCOAX, 3 = natCOAX
      std::vector<int>   GetDetectorTypes() const { return detectorTypes; }
      // get mass [g]
      std::vector<int>   GetMass()          const { return mass; }
      // get total volume, active volume, dead volume [cm^3]
      std::vector<float> GetVolume()        const;
      std::vector<float> GetActiveVolume()  const;
      std::vector<float> GetDeadVolume()    const;
      // get number of 76Ge isotopes [mol]
      std::vector<float> GetActiveN76Ge()   const;
      std::vector<float> GetDeadN76Ge()     const;

      protected:
      
      // detector's names
      std::vector<std::string> names;
      // detector's types in detector strings:
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
  
  // order vectors in the MaGe naming scheme: GELATIO -> MaGe
  template<typename T>
  void ReorderAsMaGe(std::vector<T>& v, std::string mode) {
    
    auto v_ = v;
    if (mode == "input") {
        int c = 0;
        for ( int i = 37; i <= 39 ; i++ ) { v[c] = v_[i]; c++; }
        for ( int i = 8 ; i <= 10 ; i++ ) { v[c] = v_[i]; c++; }
        for ( int i = 27; i <= 29 ; i++ ) { v[c] = v_[i]; c++; }
        v[c] = v_[36]; c++;
        for ( int i = 0 ; i <= 7  ; i++ ) { v[c] = v_[i]; c++; }
        for ( int i = 11; i <= 26 ; i++ ) { v[c] = v_[i]; c++; }
        for ( int i = 30; i <= 35 ; i++ ) { v[c] = v_[i]; c++; }
        return;
    }

    else if (mode == "output") {
        int c = 0;
        for ( int i = 37; i <= 39; i++ ) { v[c] = v_[i]; c++; }
        for ( int i = 0 ; i <= 36; i++ ) { v[c] = v_[i]; c++; }
        return;
    }

    else return;
  }
  
  // print vectors depending on the ordering scheme
  template<typename T>
  void Print(std::vector<T>& v, std::string mode) {
    
    if (mode == "GELATIO") {
        std::cout << "string1  ";
        for ( int i = 0 ; i <= 7 ; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring2  ";
        for ( int i = 8 ; i <= 10; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring3  ";
        for ( int i = 11; i <= 18; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring4  ";
        for ( int i = 19; i <= 26; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring5  ";
        for ( int i = 27; i <= 29; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring6  ";
        for ( int i = 30; i <= 36; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring7  ";
        for ( int i = 37; i <= 38; ++i ) std::cout << v[i] << ", ";
        std::cout << v[39] << std::endl;
        return;
    }
    else if (mode == "MaGeOutput") {
        std::cout << "string7  ";
        for ( int i = 0 ; i <= 2 ; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring1  ";
        for ( int i = 3 ; i <= 10; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring2  ";
        for ( int i = 11; i <= 13; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring3  ";
        for ( int i = 14; i <= 21; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring4  ";
        for ( int i = 22; i <= 29; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring5  ";
        for ( int i = 30; i <= 32; ++i ) std::cout << v[i] << ", ";
        std::cout << "\nstring6  ";
        for ( int i = 33; i <= 38; ++i ) std::cout << v[i] << ", ";
        std::cout << v[39] << std::endl;
        return;
    }
    else if (mode == "MaGeInput") {
        std::cout << "natCOAX  ";
        for ( int i = 0 ; i <= 2 ; ++i ) std::cout << v[i] << ", ";
        std::cout << "\n\nenrCOAX  ";
        for ( int i = 3 ; i <= 9 ; ++i ) std::cout << v[i] << ", ";
        std::cout << "\n\nBEGe     ";
        for ( int i = 10 ; i <= 17 ; ++i ) std::cout << v[i] << ", ";
        std::cout << "\n         ";
        for ( int i = 18 ; i <= 25 ; ++i ) std::cout << v[i] << ", ";
        std::cout << "\n         ";
        for ( int i = 26 ; i <= 33 ; ++i ) std::cout << v[i] << ", ";
        std::cout << "\n         ";
        for ( int i = 34 ; i <= 38 ; ++i ) std::cout << v[i] << ", ";
        std::cout << v[39] << std::endl;
        return;
    }
    else { std::cout << "Print failed.\n"; return; }
  }
}

#endif
