/* Fit2nbbLV.cxx
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 20/03/2017
 *
 */

#include "Fit2nbbLV.h"

#include <iostream>
#include <math.h>

#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TText.h"
#include "TF1.h"

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameter.h>

// ---------------------------------------------------------
Fit2nbbLV::Fit2nbbLV(std::string name) : BCModel(name.c_str()), kUseRange(false) {
   
    // define parameters
    /* [0] */  this->AddParameter("2nbb",                 2E02, 2.4E02);
    /* [1] */  this->AddParameter("2nbbLV",               0, 2E-03);
    /* [2] */  this->AddParameter("K42homLAr",            1E-04, 2.5E-04);
    /* [3] */  this->AddParameter("K40fibers",            0, 1E-00);
    /* [4] */  this->AddParameter("Bi212Tl208fibers",     0, 0.02);
    /* [5] */  this->AddParameter("Pb214Bi214fibers",     0, 1E-01);
    /* [6] */  this->AddParameter("alphaBEGe",            1.1E03, 1.8E03);
    /* [7] */  this->AddParameter("alphaCOAX",            2.6E03, 3.5E03);
    /* [8] */  this->AddParameter("K42nPlusBEGe",         0, 2E-04);
    /* [9] */  this->AddParameter("K42nPlusCOAX",         0, 3E-04);
    /* [10] */ this->AddParameter("K42pPlusBEGe",         0, 7E-01);
    /* [11] */ this->AddParameter("K42pPlusCOAX",         0, 1E-03);
    /* [12] */ this->AddParameter("Ac228holder",          0, 3.9E-04); // <--- Screening
    /* [13] */ this->AddParameter("Co60holder",           0, 1.6E-04); // <--- Screening
    /* [14] */ this->AddParameter("K40holder",            0, 2E-02);
    /* [15] */ this->AddParameter("Bi212Tl208holder",     0, 7E-04);
    /* [16] */ this->AddParameter("Pb214Bi214holder",     0, 30);
    /* [17] */ this->AddParameter("K40cable",             0, 1E00);
    /* [18] */ this->AddParameter("Bi212Tl208cables",     0, 5E-02);
    /* [19] */ this->AddParameter("Pb214Bi214cables",     0, 25);
    /* [20] */ this->AddParameter("K40minishroud",        0, 2E-01);
    /* [21] */ this->AddParameter("Pa234minishroud",      0, 2E-01);
    /* [22] */ this->AddParameter("Bi207minishroud",      0, 5E-03);
    /* [23] */ this->AddParameter("Bi207cables",          0, 7E-02);
    /* [24] */ this->AddParameter("Bi207holder",          0, 5E-02);
    /* [25] */ this->AddParameter("Pb214Bi214minishroud", 0, 30);
    /* [26] */ this->AddParameter("K42minishroudsurface", 0, 1E-03);
    /* [27] */ this->AddParameter("Pa234cables",          0, 3E-01);
    /* [28] */ this->AddParameter("Pa234holder",          0, 1E-02);
    /* [29] */ this->AddParameter("K42homLArAA",          0, 2E-03);

    // TODO: add new parameters here
    //// LEGEND
    //
    // [0] 2nbb
    // [1] 2nbbLV
    // [2] K42 in LAr
    // [3] K40 in fibers
    // [4] Bi212 + Tl208 in fibers
    // [5] Pb214 + Bi214 in fibers
    // [6] alphaBEGe
    // [7] alphaCOAX
    // [8] K42 nPlus BEGe
    // [9] K42 nPlus COAX
    // [10] K42 pPlus BEGe
    // [11] K42 pPlus COAX
    // [12] Ac228 holder
    // [13] Co60holder
    // [14] K40holder
    // [15] Bi212 + Tl208 in holder
    // [16] Pb214 + Bi214 in holder
    // [17] K40cables
    // [18] Bi212 + Tl208 in cables
    // [19] Pb214 + Bi214 in cables
    // [20] K40minishroud
    // [21] Pa234 minishroud
    // [22] Bi207minishroud
    // [23] Bi207cables
    // [24] Bi207holder
    // [25] Pb214 + Bi214 in minishroud
    // [26] K42 on minishroud surface
    // [27] Pa234 on cables
    // [28] Pa234 on holders
    // [29] K42 in LAr Above Array
    //

    // priors
    TF1 invflat("inverse-flat", "1/x^2"  , 0, 1);
    TF1 logflat("log-flat"    , "1/x"    , 0, 1);
    //TF1 linear ("linear"      , "(x<=[0])*(-x+[0])+(x>[0])*0" , 0, 1);

    // 2nbb
    this->SetPriorConstantAll();
    this->SetPrior(0, &invflat);

    // 2nbbLV
    //linear.SetParameter(0, 1.52E-05);
    //this->SetPriorGauss(1, 0, 1.38E-04);

    // holders
    this->SetPriorGauss(14, 4.3E-03, 0.9E-03); // K40
    //linear.SetParameter(0, 0.39E-03);          // Ac228
    //this->SetPrior(12, &linear);
    //linear.SetParameter(0, 0.16E-03);          // Co60
    //this->SetPrior(13, &linear);
}
// ---------------------------------------------------------
double Fit2nbbLV::LogLikelihood(const std::vector<double> & parameters) {

    double logprob = 0.;
    double f;

    for ( int i = downBin; i <= upBin; ++i ) {

        // ignore ROI
        if ( dbin[i] >= 2014 and dbin[i+1] <= 2064 ) continue;

        // BEGe
        f = parameters[0]*simBEGe[0][i] + parameters[0]*parameters[1]*n2n1*simBEGe[1][i]; // 2nbb & 2nbbLV
        f += parameters[2]*simBEGe[2][i];                           // K42homLAr
        f += parameters[3]*simBEGe[3][i];                           // K40fibers
        f += parameters[4]*(simBEGe[4][i] + BrTl*simBEGe[5][i]);    // Bi212 -> Tl208 fibers with 35.93% Br
        f += parameters[5]*(simBEGe[6][i] +      simBEGe[7][i]);    // Pb214 -> Bi214 fibers with 100.00% Br
        f += parameters[6]*simBEGe[8][i];                           // alpha
        f += parameters[8]*simBEGe[9][i];                           // K42nPlus
        f += parameters[10]*simBEGe[10][i];                         // K42pPlus
        f += parameters[12]*simBEGe[11][i];                         // Ac228holder
        f += parameters[13]*simBEGe[12][i];                         // Co60holder
        f += parameters[14]*simBEGe[13][i];                         // K40holder
        f += parameters[15]*(simBEGe[14][i] + BrTl*simBEGe[15][i]); // Bi212 -> Tl208 holder with 35.93% Br
        f += parameters[16]*(simBEGe[16][i] +      simBEGe[17][i]); // Pb214 -> Bi214 holder with 100.00% Br
        f += parameters[17]*simBEGe[18][i];                         // K40cables
        f += parameters[18]*(simBEGe[19][i] + BrTl*simBEGe[20][i]); // Bi212 -> Tl208 cables with 35.93% Br
        f += parameters[19]*(simBEGe[21][i] +      simBEGe[22][i]); // Pb214 -> Bi214 cables with 100.00% Br
        f += parameters[20]*simBEGe[23][i];                         // K40minishroud
        f += parameters[21]*simBEGe[24][i];                         // Pa234minishroud
        f += parameters[22]*simBEGe[25][i];                         // Bi207minishroud
        f += parameters[23]*simBEGe[26][i];                         // Bi207cables
        f += parameters[24]*simBEGe[27][i];                         // Bi207holder
        f += parameters[25]*(simBEGe[28][i] +      simBEGe[29][i]); // Pb214 -> Bi214 minishroud with 100.00% Br
        f += parameters[26]*simBEGe[30][i];                         // K42minishroudsurface
        f += parameters[27]*simBEGe[31][i];                         // Pa234cables
        f += parameters[28]*simBEGe[32][i];                         // Pa234holder
        f += parameters[29]*simBEGe[33][i];                         // K42homLArAA
        // TODO: update loglikelihood here

        //logprob += dataBEGe[i]*log(f) - f - BCMath::LogFact(dataBEGe[i]);
        logprob += BCMath::LogPoisson(dataBEGe[i], f);

        // COAX
        f = parameters[0]*simCOAX[0][i] + parameters[0]*parameters[1]*n2n1*simCOAX[1][i]; 
        f += parameters[2]*simCOAX[2][i];
        f += parameters[3]*simCOAX[3][i];
        f += parameters[4]*(simCOAX[4][i] + BrTl*simCOAX[5][i]);
        f += parameters[5]*(simCOAX[6][i] +      simCOAX[7][i]);
        f += parameters[7]*simCOAX[8][i];
        f += parameters[9]*simCOAX[9][i] + parameters[11]*simCOAX[10][i];
        f += parameters[12]*simCOAX[11][i];
        f += parameters[13]*simCOAX[12][i];
        f += parameters[14]*simCOAX[13][i];
        f += parameters[15]*(simCOAX[14][i] + BrTl*simCOAX[15][i]);
        f += parameters[16]*(simCOAX[16][i] +      simCOAX[17][i]);
        f += parameters[17]*simCOAX[18][i];
        f += parameters[18]*(simCOAX[19][i] + BrTl*simCOAX[20][i]);
        f += parameters[19]*(simCOAX[21][i] +      simCOAX[22][i]);
        f += parameters[20]*simCOAX[23][i];
        f += parameters[21]*simCOAX[24][i];
        f += parameters[22]*simCOAX[25][i];
        f += parameters[23]*simCOAX[26][i];
        f += parameters[24]*simCOAX[27][i];
        f += parameters[25]*(simCOAX[28][i] +      simCOAX[29][i]);
        f += parameters[26]*simCOAX[30][i];
        f += parameters[27]*simCOAX[31][i];
        f += parameters[28]*simCOAX[32][i];
        f += parameters[29]*simCOAX[33][i];
        // TODO: update loglikelihood here
        
        //logprob += dataCOAX[i]*log(f) - f - BCMath::LogFact(dataCOAX[i]);
        logprob += BCMath::LogPoisson(dataCOAX[i], f);
	}
    
    return logprob;
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
          + bestpar[6]*simBEGe[8][i]
          + bestpar[8]*simBEGe[9][i] + bestpar[10]*simBEGe[10][i]
          + bestpar[12]*simBEGe[11][i]
          + bestpar[13]*simBEGe[12][i]
          + bestpar[14]*simBEGe[13][i]
          + bestpar[15]*(simBEGe[14][i] + BrTl*simBEGe[15][i])
          + bestpar[16]*(simBEGe[16][i] +     simBEGe[17][i])
          + bestpar[17]*simBEGe[18][i]
          + bestpar[18]*(simBEGe[19][i] + BrTl*simBEGe[20][i])
          + bestpar[19]*(simBEGe[21][i] +     simBEGe[22][i])
          + bestpar[20]*simBEGe[23][i]
          + bestpar[21]*simBEGe[24][i]
          + bestpar[22]*simBEGe[25][i]
          + bestpar[23]*simBEGe[26][i]
          + bestpar[24]*simBEGe[27][i]
          + bestpar[25]*(simBEGe[28][i] +     simBEGe[29][i])
          + bestpar[26]*simBEGe[30][i]
          + bestpar[27]*simBEGe[31][i]
          + bestpar[28]*simBEGe[32][i]
          + bestpar[29]*simBEGe[33][i]
          ;
          // TODO: update fitting func here
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
          + bestpar[7]*simCOAX[8][i]
          + bestpar[9]*simCOAX[9][i] + bestpar[11]*simCOAX[10][i]
          + bestpar[12]*simCOAX[11][i]
          + bestpar[13]*simCOAX[12][i]
          + bestpar[14]*simCOAX[13][i]
          + bestpar[15]*(simCOAX[14][i] + BrTl*simCOAX[15][i])
          + bestpar[16]*(simCOAX[16][i] +      simCOAX[17][i])
          + bestpar[17]*simCOAX[18][i]
          + bestpar[18]*(simCOAX[19][i] + BrTl*simCOAX[20][i])
          + bestpar[19]*(simCOAX[21][i] +      simCOAX[22][i])
          + bestpar[20]*simCOAX[23][i]
          + bestpar[21]*simCOAX[24][i]
          + bestpar[22]*simCOAX[25][i]
          + bestpar[23]*simCOAX[26][i]
          + bestpar[24]*simCOAX[27][i]
          + bestpar[25]*(simCOAX[28][i] +      simCOAX[29][i])
          + bestpar[26]*simCOAX[30][i]
          + bestpar[27]*simCOAX[31][i]
          + bestpar[28]*simCOAX[32][i]
          + bestpar[29]*simCOAX[33][i]
          ;
        // TODO: update fitting func here
    }

    return totfnc;
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
// ---------------------------------------------------------
void Fit2nbbLV::WriteHistosOnFile(std::string path) {

    BCLog::OutSummary("Writing histograms on file...");
    auto results = this->GetBestFitParameters();
    int nbins = this->GetNbins();

    std::string filename = path + "outBATHist.root";
    TFile outDrawFile(filename.c_str(), "RECREATE");

    TH1D hDataBEGe("hDataBEGe", "hDataBEGe", nbins, &dbin[0]);
    TH1D hDataCOAX("hDataCOAX", "hDataBEGe", nbins, &dbin[0]);
    std::vector<TH1D> hSimBEGe;
    std::vector<TH1D> hSimCOAX;

    for ( unsigned int i = 0; i < simBEGe.size(); ++i ) {
        hSimBEGe.emplace_back(Form("h%iBEGe",i), Form("h%iBEGe",i), nbins, &dbin[0]);
        hSimCOAX.emplace_back(Form("h%iCOAX",i), Form("h%iCOAX",i), nbins, &dbin[0]);
    }

    // fill ROOT objects
    for ( int i = 0; i < nbins; ++i ) {
        hDataBEGe.SetBinContent(i+1, dataBEGe[i]);
        hDataCOAX.SetBinContent(i+1, dataCOAX[i]);

        for ( unsigned int j = 0; j < simBEGe.size(); ++j ) {
            hSimBEGe[j].SetBinContent(i+1, simBEGe[j][i]);
            hSimCOAX[j].SetBinContent(i+1, simCOAX[j][i]);
        }
    }

    // scale with best parameters and rename
    // 2nbb
    hSimBEGe[0].Scale(results[0]);
    hSimBEGe[0].SetName("h2nbbBEGe");
    hSimCOAX[0].Scale(results[0]);
    hSimCOAX[0].SetName("h2nbbCOAX");

    // 2nbbLV
    hSimBEGe[1].Scale(results[0]*results[1]*n2n1);
    hSimBEGe[1].SetName("h2nbbLVBEGe");
    hSimCOAX[1].Scale(results[0]*results[1]*n2n1);
    hSimCOAX[1].SetName("h2nbbLVCOAX");

    // K42homLAr
    hSimBEGe[2].Scale(results[2]);
    hSimBEGe[2].SetName("hK42homLArBEGe");
    hSimCOAX[2].Scale(results[2]);
    hSimCOAX[2].SetName("hK42homLArCOAX");

    // K40fibers
    hSimBEGe[3].Scale(results[3]);
    hSimBEGe[3].SetName("hK40fibersBEGe");
    hSimCOAX[3].Scale(results[3]);
    hSimCOAX[3].SetName("hK40fibersCOAX");

    // Bi212fibers
    hSimBEGe[4].Scale(results[4]);
    hSimBEGe[4].SetName("hBi212fibersBEGe");
    hSimCOAX[4].Scale(results[4]);
    hSimCOAX[4].SetName("hBi212fibersCOAX");

    // Tl208fibers
    hSimBEGe[5].Scale(results[4]*BrTl);
    hSimBEGe[5].SetName("hTl208fibersBEGe");
    hSimCOAX[5].Scale(results[4]*BrTl);
    hSimCOAX[5].SetName("hTl208fibersCOAX");

    // Pb214fibers
    hSimBEGe[6].Scale(results[5]);
    hSimBEGe[6].SetName("hPb214fibersBEGe");
    hSimCOAX[6].Scale(results[5]);
    hSimCOAX[6].SetName("hPb214fibersCOAX");

    // Bi214fibers
    hSimBEGe[7].Scale(results[5]);
    hSimBEGe[7].SetName("hBi214fibersBEGe");
    hSimCOAX[7].Scale(results[5]);
    hSimCOAX[7].SetName("hBi214fibersCOAX");

    // alphas
    hSimBEGe[8].Scale(results[6]);
    hSimBEGe[8].SetName("hAlphaBEGe");
    hSimCOAX[8].Scale(results[7]);
    hSimCOAX[8].SetName("hAlphaCOAX");

    // nPlus 
    hSimBEGe[9].Scale(results[8]);
    hSimBEGe[9].SetName("hK42nPlusBEGe");
    hSimCOAX[9].Scale(results[9]);
    hSimCOAX[9].SetName("hK42nPlusCOAX");

    // pPlus 
    hSimBEGe[10].Scale(results[10]);
    hSimBEGe[10].SetName("hK42pPlusBEGe");
    hSimCOAX[10].Scale(results[11]);
    hSimCOAX[10].SetName("hK42pPlusCOAX");

    // Ac228holder
    hSimBEGe[11].Scale(results[12]);
    hSimBEGe[11].SetName("hAc228holderBEGe");
    hSimCOAX[11].Scale(results[12]);
    hSimCOAX[11].SetName("hAc228holderCOAX");

    // Co60holder
    hSimBEGe[12].Scale(results[13]);
    hSimBEGe[12].SetName("hCo60holderBEGe");
    hSimCOAX[12].Scale(results[13]);
    hSimCOAX[12].SetName("hCo60holderCOAX");

    // Co60holder
    hSimBEGe[13].Scale(results[14]);
    hSimBEGe[13].SetName("hK40holderBEGe");
    hSimCOAX[13].Scale(results[14]);
    hSimCOAX[13].SetName("hK40holderCOAX");

    // Bi212holder
    hSimBEGe[14].Scale(results[15]);
    hSimBEGe[14].SetName("hBi212holderBEGe");
    hSimCOAX[14].Scale(results[15]);
    hSimCOAX[14].SetName("hBi212holderCOAX");

    // Tl208holder
    hSimBEGe[15].Scale(results[15]*BrTl);
    hSimBEGe[15].SetName("hTl208holderBEGe");
    hSimCOAX[15].Scale(results[15]*BrTl);
    hSimCOAX[15].SetName("hTl208holderCOAX");

    // Pb214holder
    hSimBEGe[16].Scale(results[16]);
    hSimBEGe[16].SetName("hPb214holderBEGe");
    hSimCOAX[16].Scale(results[16]);
    hSimCOAX[16].SetName("hPb214holderCOAX");

    // Bi214holder
    hSimBEGe[17].Scale(results[16]);
    hSimBEGe[17].SetName("hBi214holderBEGe");
    hSimCOAX[17].Scale(results[16]);
    hSimCOAX[17].SetName("hBi214holderCOAX");

    // K40cables
    hSimBEGe[18].Scale(results[17]);
    hSimBEGe[18].SetName("hK40cablesBEGe");
    hSimCOAX[18].Scale(results[17]);
    hSimCOAX[18].SetName("hK40cablesCOAX");

    // Bi212cables
    hSimBEGe[19].Scale(results[18]);
    hSimBEGe[19].SetName("hBi212cablesBEGe");
    hSimCOAX[19].Scale(results[18]);
    hSimCOAX[19].SetName("hBi212cablesCOAX");

    // Tl208cables
    hSimBEGe[20].Scale(results[18]*BrTl);
    hSimBEGe[20].SetName("hTl208cablesBEGe");
    hSimCOAX[20].Scale(results[18]*BrTl);
    hSimCOAX[20].SetName("hTl208cablesCOAX");

    // Pb214cables
    hSimBEGe[21].Scale(results[19]);
    hSimBEGe[21].SetName("hPb214cablesBEGe");
    hSimCOAX[21].Scale(results[19]);
    hSimCOAX[21].SetName("hPb214cablesCOAX");

    // Bi214cables
    hSimBEGe[22].Scale(results[19]);
    hSimBEGe[22].SetName("hBi214cablesBEGe");
    hSimCOAX[22].Scale(results[19]);
    hSimCOAX[22].SetName("hBi214cablesCOAX");

    // K40minishroud
    hSimBEGe[23].Scale(results[20]);
    hSimBEGe[23].SetName("hK40minishroudBEGe");
    hSimCOAX[23].Scale(results[20]);
    hSimCOAX[23].SetName("hK40minishroudCOAX");

    // Pa234minishroud
    hSimBEGe[24].Scale(results[21]);
    hSimBEGe[24].SetName("hPa234minishroudBEGe");
    hSimCOAX[24].Scale(results[21]);
    hSimCOAX[24].SetName("hPa234minishroudCOAX");

    // Bi207minishroud
    hSimBEGe[25].Scale(results[22]);
    hSimBEGe[25].SetName("hBi207minishroudBEGe");
    hSimCOAX[25].Scale(results[22]);
    hSimCOAX[25].SetName("hBi207minishroudCOAX");

    // Bi207cables
    hSimBEGe[26].Scale(results[23]);
    hSimBEGe[26].SetName("hBi207cablesBEGe");
    hSimCOAX[26].Scale(results[23]);
    hSimCOAX[26].SetName("hBi207cablesCOAX");

    // Bi207holder
    hSimBEGe[27].Scale(results[24]);
    hSimBEGe[27].SetName("hBi207holderBEGe");
    hSimCOAX[27].Scale(results[24]);
    hSimCOAX[27].SetName("hBi207holderCOAX");

    // Pb214minishroud
    hSimBEGe[28].Scale(results[25]);
    hSimBEGe[28].SetName("hPb214minishroudBEGe");
    hSimCOAX[28].Scale(results[25]);
    hSimCOAX[28].SetName("hPb214minishroudCOAX");

    // Bi214minishroud
    hSimBEGe[29].Scale(results[25]);
    hSimBEGe[29].SetName("hBi214minishroudBEGe");
    hSimCOAX[29].Scale(results[25]);
    hSimCOAX[29].SetName("hBi214minishroudCOAX");

    // K42minishroudsurface
    hSimBEGe[30].Scale(results[26]);
    hSimBEGe[30].SetName("K42minishroudsurfaceBEGe");
    hSimCOAX[30].Scale(results[26]);
    hSimCOAX[30].SetName("K42minishroudsurfaceCOAX");

    // Pa234cables
    hSimBEGe[31].Scale(results[27]);
    hSimBEGe[31].SetName("Pa234cablesBEGe");
    hSimCOAX[31].Scale(results[27]);
    hSimCOAX[31].SetName("Pa234cablesCOAX");

    // K42minishroudsurface
    hSimBEGe[32].Scale(results[28]);
    hSimBEGe[32].SetName("Pa234holderBEGe");
    hSimCOAX[32].Scale(results[28]);
    hSimCOAX[32].SetName("Pa234holderCOAX");

    // K42minishroudsurface
    hSimBEGe[33].Scale(results[29]);
    hSimBEGe[33].SetName("K42homLArAABEGe");
    hSimCOAX[33].Scale(results[29]);
    hSimCOAX[33].SetName("K42homLArAACOAX");

    // TODO: add new sources here
    // ...

    auto drawpdf = [&]( std::vector<TH1D>& v , TH1D& vd ,std::string type ) {

        TCanvas tmp(type.c_str(), type.c_str(), 2700, 700);
        TPad pad("pad", "pad", 0.0, 0.0, 1, 1);
        pad.SetMargin(0.032,0.015,0.1,0.05);
        tmp.cd();
        pad.Draw();
        pad.cd();

        vd.SetStats(false);
        vd.SetMarkerStyle(6);
        vd.GetXaxis()->SetRange(downBin,upBin);
        vd.GetXaxis()->SetTitle("energy [keV]");
        vd.GetYaxis()->SetTitle("counts");
        vd.GetYaxis()->SetNdivisions(10);
        vd.GetYaxis()->SetTickLength(0.01);
        vd.Draw("P0");

        for ( auto& h : v ) h.SetLineWidth(1);
        v[0].SetLineColor(kBlue);
        v[1].SetLineColor(kBlue+2);
        v[2].SetLineColor(kMagenta);
        // fibers
        v[3].SetLineColor(kGreen);
        v[4].SetLineColor(kGreen);
        v[5].SetLineColor(kGreen);
        v[6].SetLineColor(kGreen);
        v[7].SetLineColor(kGreen);
        // alphas
        v[8].SetLineColor(kRed);
        // contacts
        v[9].SetLineColor(kCyan);
        v[10].SetLineColor(kCyan);
        // holders
        v[11].SetLineColor(kYellow);
        v[12].SetLineColor(kYellow);
        v[13].SetLineColor(kYellow);
        v[14].SetLineColor(kYellow);
        v[15].SetLineColor(kYellow);
        v[16].SetLineColor(kYellow);
        v[17].SetLineColor(kYellow);
        // cables
        v[18].SetLineColor(kOrange);
        v[19].SetLineColor(kOrange);
        v[20].SetLineColor(kOrange);
        v[21].SetLineColor(kOrange);
        v[22].SetLineColor(kOrange);
        // minishroud
        v[23].SetLineColor(kGray);
        v[24].SetLineColor(kGray);
        v[25].SetLineColor(kGray);

        v[26].SetLineColor(kOrange);  // cables
        v[27].SetLineColor(kYellow);  // holder

        v[28].SetLineColor(kGray);  // minishroud
        v[29].SetLineColor(kGray);  // minishroud 
        v[30].SetLineColor(kBlack); // minishroudsurface
        v[31].SetLineColor(kOrange); // cables
        v[32].SetLineColor(kYellow); // holder
        v[33].SetLineColor(kMagenta+2); // LAr Above Array

        // TODO: define color

        for ( auto& h : v ) h.Draw("HISTSAME");

        // summed histogram
        std::string name = "hsum" + type;
        TH1D sum(name.c_str(), name.c_str(), nbins, &dbin[0]);
        for ( auto& h : v ) sum.Add(&h);
        sum.SetLineColor(kRed);
        sum.Draw("HISTSAME");
        sum.Write();

        TLegend leg;
        leg.SetX1NDC(0.87);
        leg.SetY1NDC(0.7);
        leg.SetX2NDC(0.98);
        leg.SetY2NDC(0.95);
        for ( auto& h : v ) leg.AddEntry(&h,h.GetName(),"l");
        leg.AddEntry(&sum,sum.GetName(),"l");
        std::string heading = "#bf{" + type + "}";
        leg.SetHeader(heading.c_str());
        leg.SetMargin(0.5);
        leg.Draw();

        //TText text(0.5, 0.7, type.c_str());
        //text.SetNDC();
        //text.Draw();

        pad.SetLogy();
        pad.SetGrid();
        //name = path + "/" + type + ".pdf";
        //tmp.SaveAs(name.c_str());
        name = path + "/" + type + ".C";
        tmp.SaveAs(name.c_str());

////////// create brazilian plot
        TCanvas brascan(type.c_str(), type.c_str(), 2700, 700);
        TPad pad1("pad", "pad", 0.0, 0.0, 1, 1);
        pad1.SetMargin(0.032,0.015,0.1,0.05);
        brascan.cd();
        pad1.Draw();
        pad1.cd();
        pad1.SetGrid();
        pad1.SetLogy();

        sum.SetStats(false);
        sum.GetXaxis()->SetRange(downBin,upBin);
        sum.GetXaxis()->SetTitle("energy [keV]");
        sum.GetYaxis()->SetTitle("counts");
        sum.GetYaxis()->SetNdivisions(10);
        sum.GetYaxis()->SetTickLength(0.01);

        TH1D resdraw(sum);
        TH1D datares(vd);

        for ( int i = 1; i <= nbins; ++i ) {
            resdraw.SetBinContent(i, 0);
            resdraw.SetBinError(i, sum.GetBinError(i));
            datares.SetBinContent(i, vd.GetBinContent(i) - sum.GetBinContent(i));
        }

        TLegend legb;
        legb.SetX1NDC(0.87);
        legb.SetY1NDC(0.1);
        legb.SetX2NDC(0.98);
        legb.SetY2NDC(0.95);
        legb.SetHeader(type.c_str());
        legb.SetMargin(0.5);

        for ( int i = 1; i <= nbins; ++i ) sum.SetBinError(i, 3*sum.GetBinError(i));
        sum.SetMarkerStyle(0);
        sum.SetFillColor(kOrange-3);
        sum.SetFillStyle(1001);
        sum.DrawCopy("E2");
        legb.AddEntry(sum.Clone(), "3#sigma", "f");

        for ( int i = 1; i <= nbins; ++i ) sum.SetBinError(i, 2*sum.GetBinError(i)/3);
        sum.SetFillColor(kYellow);
        sum.DrawCopy("E2 SAME");
        legb.AddEntry(sum.Clone(), "2#sigma", "f");

        for ( int i = 1; i <= nbins; ++i ) sum.SetBinError(i, 0.5*sum.GetBinError(i));
        sum.SetFillColor(kGreen);
        sum.DrawCopy("E2 SAME");
        legb.AddEntry(sum.Clone(), "1#sigma", "f");

        for ( int i = 1; i <= nbins; ++i ) sum.SetBinError(i, 0.00001);
        sum.SetFillStyle(0);
        sum.SetFillColor(0);
        sum.DrawCopy("E SAME");
        legb.AddEntry(sum.Clone(), "Model", "l");

        vd.DrawCopy("P0 SAME");
        legb.AddEntry(sum.Clone(), "Data", "p");
        //text.Draw();
        legb.Draw();

        name = path + "/brasilian" + type + ".C";
        brascan.SaveAs(name.c_str());

///////// create residuals plot
        for ( int i = 1; i <= nbins; ++i ) resdraw.SetBinError(i, 3*resdraw.GetBinError(i));
        resdraw.SetMarkerStyle(0);
        resdraw.SetFillColor(kOrange-3);
        resdraw.SetFillStyle(1001);
        pad1.Draw();
        pad1.cd();
        pad1.SetLogy(false);
        resdraw.DrawCopy("E2");

        for ( int i = 1; i <= nbins; ++i ) resdraw.SetBinError(i, 2*resdraw.GetBinError(i)/3);
        resdraw.SetFillColor(kYellow);
        resdraw.DrawCopy("E2 SAME");

        for ( int i = 1; i <= nbins; ++i ) resdraw.SetBinError(i, 0.5*resdraw.GetBinError(i));
        resdraw.SetFillColor(kGreen);
        resdraw.DrawCopy("E2 SAME");

        for ( int i = 1; i <= nbins; ++i ) resdraw.SetBinError(i, 0.00001);
        resdraw.SetFillStyle(0);
        resdraw.SetFillColor(0);
        resdraw.DrawCopy("E SAME");

        datares.DrawCopy("P0 SAME");
        //text.Draw();
        legb.Draw();

        name = path + "/residuals" + type + ".C";
        brascan.SaveAs(name.c_str());
    };

    drawpdf( hSimBEGe, hDataBEGe, "BEGe" );
    drawpdf( hSimCOAX, hDataCOAX, "COAX" );

    hDataBEGe.Write();
    hDataCOAX.Write();
    for ( unsigned int i = 0; i < hSimBEGe.size(); ++i ) {
        hSimBEGe[i].Write();
        hSimCOAX[i].Write();
    }

    outDrawFile.Close();

    return;
}
