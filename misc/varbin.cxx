#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH1D.h"

int main() {
    
    std::string path = std::string(std::getenv("GERDACPTDIR")) + "/out/sumData.root";
    TFile file(path.c_str(), "READ");

    TH1D* h = dynamic_cast<TH1D*>(file.Get("energyBEGeAll"));

    std::vector<int> avoid = { 568, 572, 580, 584, 608, 612, 908, 
                               912, 968, 972, 1000, 1004, 1060, 1064, 
                               1120, 1172, 1176, 1236, 1240, 1332, 1460, 
                               1464, 1524, 1528, 1764, 2204, 2612, 2616 };

    double dbin[1848];
    int k = 0, i = 0;
    while (1) { 
        if ( std::find( avoid.begin(), avoid.end(), k) == avoid.end() ) { 
            dbin[i] = k; 
            i++; 
        }
        k += 4;
        if ( k > 7500 ) break;
    }

    //for( int i = 0; i < 1848; i++ ) std::cout << dbin[i] << ' ';
    
    auto hnew = h->Rebin(1827, "hnew", dbin);

    TFile outfile ("rebinned.root","RECREATE");
    hnew->Write();

    return 0;
}
