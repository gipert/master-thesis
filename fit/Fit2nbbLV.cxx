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
   
    // define parameters
    this->AddParameter("2nbb", 0, 500);
    this->AddParameter("2nbbLV", 0, 1);

    //fake
    this->AddParameter("fake1", 0, 1);
    this->AddParameter("fake2", 0, 1);

    // priors
    this->SetPriorGauss(0,217,11);
    this->SetPriorConstant(1);
    this->SetPriorConstant(2);    
    this->SetPriorConstant(3);    
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
    //int last = parameters.size()-1;
    
    for ( int i = downBin; i < upBin; ++i ) {

        // BEGe
        f = parameters[0]*simBEGe[0][i] + parameters[0]*parameters[1]*n2n1*simBEGe[1][i]; 
        //for ( int j = 2; j < last; ++j ) f += parameters[j]*simBEGe[j][i];
        logprob += dataBEGe[i]*log(f) - f;

        // COAX
        f = parameters[0]*simCOAX[0][i] + parameters[0]*parameters[1]*n2n1*simCOAX[1][i]; 
        //for ( int j = 2; j < last; ++j ) f += parameters[j]*simCOAX[j][i];
        logprob += dataCOAX[i]*log(f) - f;
	}

    return logprob;
}
