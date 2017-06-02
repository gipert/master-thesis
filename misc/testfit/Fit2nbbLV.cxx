/* Fit2nbbLV.cxx
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 20/03/2017
 *
 */

#include "Fit2nbbLV.h"

#include <iostream>
#include <math.h>

#include "TF1.h"

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameter.h>

// ---------------------------------------------------------
Fit2nbbLV::Fit2nbbLV(std::string name) : BCModel(name.c_str()), kUseRange(false) {

    // define parameters
    /* [0] */  this->AddParameter("2nbb", 0, 2);
    /* [1] */  this->AddParameter("2nbbLV", 0, 0.00001);
    // priors
    TF1 invflat("inverse-flat", "1/x^2"  , 0, 1);
    //TF1 logflat("log-flat"    , "1/x"    , 0, 1);
    //TF1 linear ("linear"      , "(x<=[0])*(-x+[0])+(x>[0])*0" , 0, 1);

    // 2nbb
    this->SetPriorConstantAll();
    //this->SetPrior(0, &invflat);

    // 2nbbLV
    //linear.SetParameter(0, 1.52E-05);
    //this->SetPriorGauss(1, 0, 1.38E-04);
}
// ---------------------------------------------------------
double Fit2nbbLV::LogLikelihood(const std::vector<double> & parameters) {

    double logprob = 0.;
    double f;

    for ( int i = downBin; i <= upBin; ++i ) {

        // BEGe
        f += parameters[0]*simBEGe[0][i] + parameters[0]*parameters[1]*n2n1*simBEGe[1][i]; // 2nbb & 2nbbLV
        logprob += BCMath::LogPoisson(dataBEGe[i], f);
    }
    return logprob;
}
// ---------------------------------------------------------
void Fit2nbbLV::SetBinning(std::vector<double>& v) {
    dbin = v;
    downBin = 0;
    upBin = v.size()-2; // this because the last entry is fake,
                        // it contains only the upper bound of 
                        // the last bin
                        // => nBins = ubin.size()-1-1
    return;
}
// ---------------------------------------------------------
void Fit2nbbLV::SetFitRange(double down, double up) {

    if (dbin.empty()) {
        std::cerr << "Error : you must call Fit2nbbLV::SetBinning first!\n" << std::flush;
        return;
    }

    int idown = 0; 
    int iup   = 0;
    int dbinsize = dbin.size();

    // maybe this needs some pre-testing?
    for ( int i = 0; i < dbinsize; ++i ) {
        if ( (double)dbin[i] <  down ) idown++;
        if ( (double)dbin[i] <= up   ) iup++;
    }

    if ( (idown != 0 or iup != 0) and iup > idown ) {
        kUseRange = true;
        downBin = idown;
        upBin = iup;

        std::cout << "Summary : Using range from dbin[" << idown << "] = " << dbin[idown]
                  << "keV to dbin[" << iup << "] = " << dbin[iup] << "keV\n";
    }

    else std::cout << "Error : Fit2nbbLV::SetFitRange: something went wrong...\n";

    return;
}
