/* Fit2nbbLV.h
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

#include <BAT/BCMath.h>
#include <BAT/BCLog.h>

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
        hDataCOAX.SetBinContent(i+1, dataBEGe[i]);
        
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
    
    // K42
    hSimBEGe[2].Scale(results[2]);
    hSimBEGe[2].SetName("hK42homLArBEGe");
    hSimCOAX[2].Scale(results[2]);
    hSimCOAX[2].SetName("hK42homLArCOAX");
    
    // K40
    hSimBEGe[3].Scale(results[3]);
    hSimBEGe[3].SetName("hK40onFiberShroudBEGe");
    hSimCOAX[3].Scale(results[3]);
    hSimCOAX[3].SetName("hK40onFiberShroudCOAX");
     
    // Bi212
    hSimBEGe[4].Scale(results[4]);
    hSimBEGe[4].SetName("hBi212onFiberShroudBEGe");
    hSimCOAX[4].Scale(results[4]);
    hSimCOAX[4].SetName("hBi212onFiberShroudCOAX");
    
    // Tl208
    hSimBEGe[5].Scale(results[4]*BrTl);
    hSimBEGe[5].SetName("hTl208onFiberShroudBEGe");    
    hSimCOAX[5].Scale(results[4]*BrTl);
    hSimCOAX[5].SetName("hTl208onFiberShroudCOAX");
    
    // Pb214
    hSimBEGe[6].Scale(results[5]);
    hSimBEGe[6].SetName("hPb214onFiberShroudBEGe");
    hSimCOAX[6].Scale(results[5]);
    hSimCOAX[6].SetName("hPb214onFiberShroudCOAX");
    
    // Bi214
    hSimBEGe[7].Scale(results[5]);
    hSimBEGe[7].SetName("hBi214onFiberShroudBEGe");
    hSimCOAX[7].Scale(results[5]);
    hSimCOAX[7].SetName("hBi214onFiberShroudCOAX");

    // alphas
    hSimBEGe[8].Scale(results[6]);
    hSimBEGe[8].SetName("hAlphaBEGe");
    hSimCOAX[8].Scale(results[7]);
    hSimCOAX[8].SetName("hAlphaCOAX");
    
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
        vd.Draw("P");
        
        for ( auto& h : v ) h.SetLineWidth(1);
        v[0].SetLineColor(kBlue);
        v[1].SetLineColor(kBlue+2);  
        v[2].SetLineColor(kGreen);
        v[3].SetLineColor(kGreen+1);
        v[4].SetLineColor(kGreen+2);
        v[5].SetLineColor(kGreen+3);
        v[6].SetLineColor(kGreen+4);
        v[7].SetLineColor(kBlack);
        v[8].SetLineColor(kRed+2);

        for ( auto& h : v ) h.Draw("HISTSAME");
        
        // summed histogram
        std::string name = "hsum" + type;
        TH1D sum(name.c_str(), name.c_str(), nbins, &dbin[0]);
        for ( auto& h : v ) sum.Add(&h);
        sum.SetLineColor(kRed);
        sum.Draw("HISTSAME");
        sum.Write();

        TLegend leg(0.87,0.64,0.98,0.95);
        for ( auto& h : v ) leg.AddEntry(&h,h.GetName(),"l");
        leg.AddEntry(&sum,sum.GetName(),"l");
        leg.Draw();
        
        pad.SetLogy();
        pad.SetGrid();
        //name = path + "/" + type + ".pdf";
        //tmp.SaveAs(name.c_str());
        name = path + "/" + type + ".C";
        tmp.SaveAs(name.c_str());

        // create brazilian plot
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

        for ( int i = 1; i <= nbins; ++i ) sum.SetBinError(i, 2*sum.GetBinError(i));
        sum.SetMarkerStyle(0);
        sum.SetFillColor(kYellow);
        sum.SetFillStyle(3001);
        sum.DrawCopy("E2");

        for ( int i = 1; i <= nbins; ++i ) sum.SetBinError(i, 0.5*sum.GetBinError(i));
        sum.SetFillColor(kGreen);
        sum.DrawCopy("E2 SAME");
        
        for ( int i = 1; i <= nbins; ++i ) sum.SetBinError(i, 0.0001);
        sum.SetFillStyle(0);
        sum.SetFillColor(0);
        sum.DrawCopy("E SAME");

        vd.DrawCopy("P SAME");
        
        name = path + "/brasilian" + type + ".C";
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
