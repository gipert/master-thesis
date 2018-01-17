#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TApplication.h"

void AddResolutionBEGe( TH1 * h );

int main( int argc, char** argv) {

    TFile f("../out/processed/K42homLAr.root");
    TH1::AddDirectory(false);
    TH1* h = (TH1*)f.Get("h0");
    AddResolutionBEGe(h);
    TApplication app("app",&argc,argv);

    TFile out("outfile.root","recreate");
    h->Write();
    return 0;
}
