#include <iostream>
#include <vector>
#include <string>
#include "DetectorSet.h"
#include "DataReader.h"

int main() {

    std::string rootpath = std::string(std::getenv("GERDACPTDIR"));
    GERDA::DetectorSet set("GELATIO");

    int BEGeMass = 0, EnrCoaxMass = 0, NatCoaxMass = 0;

    for ( int i = 0; i < 40; i++ ) {
        if ( set.GetDetectorTypes()[i] == 1 ) BEGeMass += set.GetMass()[i];
        if ( set.GetDetectorTypes()[i] == 2 ) EnrCoaxMass += set.GetMass()[i];
        if ( set.GetDetectorTypes()[i] == 3 ) NatCoaxMass += set.GetMass()[i];
    }

    std::cout << "The total mass of the BEGe detectors is: "    << BEGeMass/1E03 << '\n';
    std::cout << "The total mass of the EnrCoax detectors is: " << EnrCoaxMass/1E03 << '\n';
    std::cout << "The total mass of the NatCoax detectors is: " << NatCoaxMass/1E03 << '\n';

    GERDA::DataReader reader( rootpath + "/misc/paths.txt", false, "GELATIO");

    // get detectorStatusMap for selected runs
    auto dsm = reader.GetDetectorStatusMap();

    // get live times saved in out/results.dat and store [s]
    std::string path = rootpath + "/data/sumData.dat";
    std::ifstream timeFile(path.c_str());
    std::map<int,int> timeMap;
    int runID, time;
    while ( timeFile >> runID >> time ) timeMap.insert(std::make_pair(runID,time));

    ////// calculate the total time each detector is ON --> detector status = 0
    // totalTime follows the GELATIO scheme
    std::vector<int> totalTime(40, 0);
    for ( const auto& i : dsm ) {
        for ( int j = 0; j < 40; ++j ) {
            if ( i.second[j] == 0 ) totalTime[j] += timeMap[i.first];
        }
    }

    double BEGeExp = 0, EnrCoaxExp = 0, NatCoaxExp = 0;
    for ( int i = 0; i < 40; i++ ) {
        if ( set.GetDetectorTypes()[i] == 1 ) BEGeExp    += (set.GetMass()[i]/1E3)*(totalTime[i]/31536000);
        if ( set.GetDetectorTypes()[i] == 2 ) EnrCoaxExp += (set.GetMass()[i]/1E3)*(totalTime[i]/31536000);
        if ( set.GetDetectorTypes()[i] == 3 ) NatCoaxExp += (set.GetMass()[i]/1E3)*(totalTime[i]/31536000);
    }

    std::cout << "The total exposure of the BEGe detectors is: "    << BEGeExp << '\n';
    std::cout << "The total exposure of the EnrCoax detectors is: " << EnrCoaxExp << '\n';
    std::cout << "The total exposure of the NatCoax detectors is: " << NatCoaxExp << '\n';

    return 0;
}
