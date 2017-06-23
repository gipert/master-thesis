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

    std::cout << "The total mass of the BEGe detectors is: "    << BEGeMass/1E03 << " kg\n";
    std::cout << "The total mass of the EnrCoax detectors is: " << EnrCoaxMass/1E03 << " kg\n";
    std::cout << "The total mass of the NatCoax detectors is: " << NatCoaxMass/1E03 << " kg\n";

    GERDA::DataReader reader( rootpath + "/misc/paths.txt", false, "GELATIO");

    // get detectorStatusMap for selected runs
    auto dsm = reader.GetDetectorStatusMap();

    // get live times saved in out/results.dat and store [s]
    std::string path = rootpath + "/data/sumData.dat";
    std::ifstream timeFile(path.c_str());
    std::map<int,int> timeMap;
    int runID, time;
    while ( timeFile >> runID >> time ) timeMap.insert(std::make_pair(runID,time));

    // calculate total livetime
    int totliv = 0;
    for ( const auto& i : timeMap ) totliv += i.second;
    std::cout << "The total live time is " << (double)totliv/31536000 << " yr\n" << std::flush;

    ////// calculate the total time each detector is ON --> detector status = 0
    // totalTime follows the GELATIO scheme
    std::vector<double> totalTime(40, 0);
    for ( const auto& i : dsm ) { // <--- loop on runIDs
        for ( int j = 0; j < 40; ++j ) { // <--- loop in detectors
            if ( i.second[j] == 0 ) totalTime[j] += timeMap[i.first];
        }
    }
    // calculate livetimes
    double ltBEGe = 0, ltEnrCOAX = 0, ltNatCOAX = 0;
    for ( int i = 0; i < 40; i++ ) {
        if      ( set.GetDetectorTypes()[i] == 1 ) ltBEGe += totalTime[i];
        else if ( set.GetDetectorTypes()[i] == 2 ) ltEnrCOAX += totalTime[i];
        else if ( set.GetDetectorTypes()[i] == 3 ) ltNatCOAX += totalTime[i];
    }

    std::cout << "BEGe livetime: " << ltBEGe/31536000 << " yr\n";
    std::cout << "EnrCOAX livetime: " << ltEnrCOAX/31536000 << " yr\n";
    std::cout << "NatCOAX livetime: " << ltNatCOAX/31536000 << " yr\n";

    double BEGeExp = 0, EnrCoaxExp = 0, NatCoaxExp = 0;
    for ( int i = 0; i < 40; i++ ) {
        if ( set.GetDetectorTypes()[i] == 1 ) BEGeExp    += ((double)set.GetMass()[i]/1E3)*(totalTime[i]/31536000);
        if ( set.GetDetectorTypes()[i] == 2 ) EnrCoaxExp += ((double)set.GetMass()[i]/1E3)*(totalTime[i]/31536000);
        if ( set.GetDetectorTypes()[i] == 3 ) NatCoaxExp += ((double)set.GetMass()[i]/1E3)*(totalTime[i]/31536000);
    }

    std::cout << "The total exposure of the BEGe detectors is: "    << BEGeExp << " kg•yr\n";
    std::cout << "The total exposure of the EnrCoax detectors is: " << EnrCoaxExp << " kg•yr\n";
    std::cout << "The total exposure of the NatCoax detectors is: " << NatCoaxExp << " kg•yr\n";

    return 0;
}
