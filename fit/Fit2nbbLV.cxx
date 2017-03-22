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
   
    this->AddParameter("2nbb", 0, 0.01, "2nbb activity [Bq/kg]");
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
    double p2nbb = parameters[0];
    
    for ( int i = downBin; i < upBin; ++i ) {
        // BEGe
        logprob += dataBEGe[i]*log(p2nbb*simBEGe[0][i]) - p2nbb*simBEGe[0][i] - BCMath::LogFact(dataBEGe[i]);
        // COAX
        logprob += dataCOAX[i]*log(p2nbb*simCOAX[0][i]) - p2nbb*simCOAX[0][i] - BCMath::LogFact(dataCOAX[i]);
	}

    return logprob;
}

// ---------------------------------------------------------
// double Fit2nbbLV::LogAPrioriProbability(const std::vector<double> & parameters) {
// 	// This method returns the logarithm of the prior probability for the
// 	// parameters p(parameters).

// 	// You need not overload this function, if you are using built-in
// 	// priors through the function SetPriorGauss, SetPriorConstant, etc.
// }
// ---------------------------------------------------------

