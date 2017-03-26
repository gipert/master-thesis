/* Fit2nbbLV.h
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 20/03/2017
 *
 */

#include "Fit2nbbLV.h"

#include <iostream>
#include <math.h>

#include <BAT/BCMath.h>

// ---------------------------------------------------------
Fit2nbbLV::Fit2nbbLV(std::string name) : BCModel(name.c_str()), kUseRange(false) {
   
    this->AddParameter("2nbb", 0, 2000, "log2*N_A/Tbb");
    this->AddParameter("2nbbLV", 0, 2000, "log2*N_A/TbbLV");
    this->AddParameter("2nbb1", 0, 2000, "log2*N_A/Tbb");
    this->AddParameter("2nbb2", 0, 2000, "log2*N_A/Tbb");
    this->SetPriorConstantAll();    
}
// ---------------------------------------------------------
void Fit2nbbLV::SetBinning(std::vector<int>& v) {
    ubin = v;
    downBin = 0;
    upBin = v.size()-1;
    return;
}
// ---------------------------------------------------------
void Fit2nbbLV::SetFitRange(double down, double up) {
    
    if (ubin.empty()) {
        std::cerr << "Error: you must call Fit2nbbLV::SetBinning first!\n" << std::flush;
        return;
    }

    int idown = 0; 
    int iup   = 0;
    int ubinsize = ubin.size();
    
    // maybe this needs some pre-testing
    for ( int i = 0; i < ubinsize; ++i ) {
        if ( (double)ubin[i] <= down ) idown++;
        if ( (double)ubin[i] <  up   ) iup++;
    }

    if ( idown != 0 and iup != 0 and iup > idown ) {
        kUseRange = true;
        downBin = idown;
        upBin = iup;
    }

    else std::cout << "Fit2nbbLV::SetFitRange: something went wrong...\n";

    return;
}
// ---------------------------------------------------------
double Fit2nbbLV::LogLikelihood(const std::vector<double> & parameters) {
    
    double logprob = 0.;
    double f;
    int size = 2; // <------------ da sistemare poi
    
    for ( int i = downBin; i < upBin; ++i ) {

        // BEGe
        f = 0; for ( int j = 0; j < size; ++j ) f += parameters[j]*simBEGe[j][i];
        logprob += dataBEGe[i]*log(f) - f;

        // COAX
        f = 0; for ( int j = 0; j < size; ++j ) f += parameters[j]*simCOAX[j][i];
        logprob += dataCOAX[i]*log(f) - f;
	}

    return logprob;
}
