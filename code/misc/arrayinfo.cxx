#include <iostream>

#include "DetectorSet.h"

int main() {

    GERDA::DetectorSet set("MaGeInput");

        std::cout << "--------------------------------\n";
        std::cout << "Name" << '\t' << "Type" << '\t' << "Mass" << '\t' << "AMass\n";
        std::cout << "--------------------------------\n";
        for ( int i = 0; i < 40; ++i ) {
            std::cout << set.GetDetectorNames()[i] << '\t';

            if ( set.GetDetectorTypes()[i] == 1 ) std::cout << "BEGe";
            if ( set.GetDetectorTypes()[i] == 2 ) std::cout << "EnrCoax";
            if ( set.GetDetectorTypes()[i] == 3 ) std::cout << "NatCoax";
            std::cout << '\t';

            std::cout << set.GetMass()[i] << '\t';

            std::cout << set.GetActiveMass()[i] << '\t';
            std::cout << std::endl;

        }
        std::cout << "--------------------------------\n";

    return 0;
}
