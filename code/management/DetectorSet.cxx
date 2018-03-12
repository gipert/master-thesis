/* DetectorSet.cxx
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

#include "DetectorSet.h"

GERDA::DetectorSet::DetectorSet( std::string naming ) :

    names { /* SECRET */ },

    detectorTypes { /* SECRET */ }, // 1 BEGe, 2 enrCoax, 3 natCoax

    mass { /* SECRET */ }, // g

    fractionAV { /* SECRET */ }, // if 0.000 it's N.A.

    fractionEnr { /* SECRET */ },

    meanMolarMass { /* SECRET */ }
{
    std::string opt;
    if      (naming == "MaGeInput")  opt = "input";
    else if (naming == "MaGeOutput") opt = "output";
    else if (naming == "GELATIO")    opt = "gelatio";
    else std::cout << naming << ": unknown ordering scheme.\n";

        GERDA::ReorderAsMaGe(names, opt);
        GERDA::ReorderAsMaGe(detectorTypes, opt);
        GERDA::ReorderAsMaGe(mass, opt);
        GERDA::ReorderAsMaGe(fractionAV, opt);
        GERDA::ReorderAsMaGe(fractionEnr, opt);
        GERDA::ReorderAsMaGe(meanMolarMass, opt);
}

// -------------------------------------------------------------------------------
std::vector<float> GERDA::DetectorSet::GetVolume() const {

    std::vector<float> volume;
    for ( int i = 0; i < 40; ++i ) {
        if ( detectorTypes[i] == 1 or 
             detectorTypes[i] == 2    ) volume.push_back(((float)mass.at(i))/enrGeDensity);
        else                            volume.push_back(((float)mass.at(i))/natGeDensity);
    }
    return volume;
}
// -------------------------------------------------------------------------------
std::vector<float> GERDA::DetectorSet::GetActiveMass() const {

    std::vector<float> volume;
    for ( int i = 0; i < 40; ++i ) {
        if ( detectorTypes[i] == 1 or
             detectorTypes[i] == 2    ) volume.push_back(((float)mass.at(i))*fractionAV.at(i));
        else                            volume.push_back(((float)mass.at(i))*fractionAV.at(i));
    }
    return volume;
}
// -------------------------------------------------------------------------------
std::vector<float> GERDA::DetectorSet::GetActiveVolume() const {

    std::vector<float> volume;
    for ( int i = 0; i < 40; ++i ) {
        if ( detectorTypes[i] == 1 or
             detectorTypes[i] == 2    ) volume.push_back(((float)mass.at(i))*fractionAV.at(i)/enrGeDensity);
        else                            volume.push_back(((float)mass.at(i))*fractionAV.at(i)/natGeDensity);
    }
    return volume;
}
// -------------------------------------------------------------------------------
std::vector<float> GERDA::DetectorSet::GetDeadVolume() const {

    std::vector<float> volume;
    for ( int i = 0; i < 40; ++i ) {
        if ( detectorTypes[i] == 1 or
             detectorTypes[i] == 2    ) volume.push_back(((float)mass.at(i))*(1-fractionAV.at(i))/enrGeDensity);
        else                            volume.push_back(((float)mass.at(i))*(1-fractionAV.at(i))/natGeDensity);
    }
    return volume;
}
// -------------------------------------------------------------------------------
std::vector<float> GERDA::DetectorSet::GetActiveN76Ge() const {

    std::vector<float> N;
    std::vector<float> V = this->GetVolume();
    for ( int i = 0; i < 40; ++i ) {
        if ( detectorTypes[i] == 1 or
             detectorTypes[i] == 2    ) N.push_back(enrGeDensity*fractionAV[i]*fractionEnr[i]*V[i]/meanMolarMass[i]);
        else                            N.push_back(natGeDensity*fractionAV[i]*fractionEnr[i]*V[i]/meanMolarMass[i]);
    }
    return N;
}
// -------------------------------------------------------------------------------
std::vector<float> GERDA::DetectorSet::GetDeadN76Ge() const {

    std::vector<float> N;
    std::vector<float> V = this->GetVolume();
    for ( int i = 0; i < 40; ++i ) {
        if ( detectorTypes[i] == 1 or
             detectorTypes[i] == 2    ) N.push_back(enrGeDensity*(1-fractionAV[i])*fractionEnr[i]*V[i]/meanMolarMass[i]);
        else                            N.push_back(natGeDensity*(1-fractionAV[i])*fractionEnr[i]*V[i]/meanMolarMass[i]);
    }
    return N;
}
