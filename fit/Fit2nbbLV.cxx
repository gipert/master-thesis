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
    /* [0] */ this->AddParameter("2nbb", 200, 280);
    /* [1] */ this->AddParameter("2nbbLV", 0, 0.001);
    /* [2] */ this->AddParameter("K42homLAr", 0, 0.0003);
    /* [3] */ this->AddParameter("K40onFiberShroud", 0, 2);
    /* [4] */ this->AddParameter("Bi212onFiberShroud", 0, 0.1);
    /* [5] */ this->AddParameter("Bi214onFiberShroud", 0, 0.1);
    /* [6] */ this->AddParameter("alphaBEGe", 0, 2000);
    /* [7] */ this->AddParameter("alphaCOAX", 0, 2600);
    //
    //// LEGEND
    //
    // [0] 2nbb
    // [1] 2nbbLV
    // [2] K42 in LAr
    // [3] K40 in fibers
    // [4] Pb212 + Tl208 in fibers
    // [5] Pb214 + Bi214 in fibers
    // [6] alphaBEGe
    // [7] alphaCOAX
    //

    // priors
    /*this->SetPriorGauss(0,217,11);
    this->SetPriorConstant(1);
    this->SetPriorConstant(2);    
    this->SetPriorConstant(3);    
    this->SetPriorConstant(4);    
    this->SetPriorConstant(5);
    this->SetPriorConstant(6);
    this->SetPriorConstant(7);*/
    this->SetPriorConstantAll();
}
// ---------------------------------------------------------
std::vector<double> Fit2nbbLV::GetFittedFncBEGe(std::vector<double>& bestpar) {
    
    int nbins = this->GetNbins();
    std::vector<double> totfnc(nbins, 0);

    for ( int i = downBin; i <= upBin; ++i ) {

        totfnc[i] = 
            // BEGe
            bestpar[0]*simBEGe[0][i] + bestpar[0]*bestpar[1]*n2n1*simBEGe[1][i]
          + bestpar[2]*simBEGe[2][i] + bestpar[3]*simBEGe[3][i]
          + bestpar[4]*(simBEGe[4][i] + BrTl*simBEGe[5][i])
          + bestpar[5]*(simBEGe[6][i] +      simBEGe[7][i])
          + bestpar[6]*simBEGe[8][i];
    }

    return totfnc;
}
// ---------------------------------------------------------
std::vector<double> Fit2nbbLV::GetFittedFncCOAX(std::vector<double>& bestpar) {
    
    int nbins = this->GetNbins();
    std::vector<double> totfnc(nbins, 0);

    for ( int i = downBin; i <= upBin; ++i ) {

        totfnc[i] = 
            // COAX
            bestpar[0]*simCOAX[0][i] + bestpar[0]*bestpar[1]*n2n1*simCOAX[1][i]
          + bestpar[2]*simCOAX[2][i] + bestpar[3]*simCOAX[3][i]
          + bestpar[4]*(simCOAX[4][i] + BrTl*simCOAX[5][i])
          + bestpar[5]*(simCOAX[6][i] +      simCOAX[7][i])
          + bestpar[7]*simCOAX[8][i];
    }

    return totfnc;
}
// ---------------------------------------------------------
void Fit2nbbLV::SetBinning(std::vector<double>& v) {
    dbin = v;
    downBin = 0;
    upBin = v.size()-2;
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
    
    // maybe this needs some pre-testing
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
// ---------------------------------------------------------
double Fit2nbbLV::LogLikelihood(const std::vector<double> & parameters) {
    
    double logprob = 0.;
    double f;
    
    for ( int i = downBin; i <= upBin; ++i ) {

        // BEGe
        f = parameters[0]*simBEGe[0][i] + parameters[0]*parameters[1]*n2n1*simBEGe[1][i]; 
        for ( int j = 2; j <= 3; ++j ) f += parameters[j]*simBEGe[j][i];
        f += parameters[4]*(simBEGe[4][i] + BrTl*simBEGe[5][i]);
        f += parameters[5]*(simBEGe[6][i] +      simBEGe[7][i]);
        f += parameters[6]*simBEGe[8][i];
        
        //logprob += dataBEGe[i]*log(f) - f - BCMath::LogFact(dataBEGe[i]);
        logprob += BCMath::LogPoisson(dataBEGe[i], f);

        // COAX
        f = parameters[0]*simCOAX[0][i] + parameters[0]*parameters[1]*n2n1*simCOAX[1][i]; 
        for ( int j = 2; j <= 3; ++j ) f += parameters[j]*simCOAX[j][i];
        f += parameters[4]*(simCOAX[4][i] + BrTl*simCOAX[5][i]);
        f += parameters[5]*(simCOAX[6][i] +      simCOAX[7][i]);
        f += parameters[7]*simCOAX[8][i];
        
        //logprob += dataCOAX[i]*log(f) - f - BCMath::LogFact(dataCOAX[i]);
        logprob += BCMath::LogPoisson(dataCOAX[i], f);
	}
    
    return logprob;
}
