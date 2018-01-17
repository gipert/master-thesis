#include <iostream>
#include <string>

#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TRandom3.h"

// h must be an histogram with 'true' entries
void AddResolutionBEGe( TH1 * h ) {

    // retrieve resolution curves
    std::string path(std::getenv("GERDACPTDIR"));
    path += "/misc/res_curves.root";
    TFile inFile(path.c_str(),"READ");
    TF1 * resBEGe = dynamic_cast<TF1*>(inFile.Get("ResCurveBege"));
    if ( !resBEGe ) { std::cout << "Zombie BEGe resolution curve!\n"; return; }

    TH1 * h_tmp = dynamic_cast<TH1*>(h->Clone());
    h->Reset();
    TRandom3 gen(0);
    int binCont;
    double binCenter;
    double u;
    //bool warn = false;
    for ( int i = 1; i < h_tmp->GetNbinsX(); ++i ) {
        binCont = h_tmp->GetBinContent(i);
        //if ( binCont != (int)binCont and warn == false ) {
        //    warn = true;
        //    std::cout << "Warning: found non-integer content in bin " << i << "!\n";
        //}
        binCenter = h_tmp->GetBinCenter(i);
        for ( int j = 0; j < binCont; ++j ) {
            u = gen.Gaus(0, resBEGe->Eval(binCenter)/2.355); // FWHM -> sigma
            h->Fill(binCenter + u);
        }
    }

    delete h_tmp;
    return;
}

void AddResolutionCOAX( TH1 * h ) {

    // retrieve resolution curves
    std::string path(std::getenv("GERDACPTDIR"));
    path += "/misc/res_curves.root";
    TFile inFile(path.c_str(),"READ");
    TF1 * resCoax = dynamic_cast<TF1*>(inFile.Get("ResCurveCoax"));
    if ( !resCoax ) { std::cout << "Zombie BEGe resolution curve!\n"; return; }

    TH1 * h_tmp = dynamic_cast<TH1*>(h->Clone());
    h->Reset();
    TRandom3 gen(0);
    int binCont;
    double binCenter;
    double u;
    //bool warn = false;
    for ( int i = 1; i < h_tmp->GetNbinsX(); ++i ) {
        binCont = h_tmp->GetBinContent(i);
        //if ( binCont != (int)binCont and warn == false ) {
        //    warn = true;
        //    std::cout << "Warning: found non-integer content in bin " << i << "!\n";
        //}
        binCenter = h_tmp->GetBinCenter(i);
        for ( int j = 0; j < binCont; ++j ) {
            u = gen.Gaus(0, resCoax->Eval(binCenter)/2.355); // FWHM -> sigma
            h->Fill(binCenter + u);
        }
    }

    delete h_tmp;
    return;
}
