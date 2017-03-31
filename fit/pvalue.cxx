#include <iostream>
#include <vector>
#include <random>

#include "ProgressBar.h"

#include "BAT/BCMath.h"
#include "Fit2nbbLV.h"

double GetPValue(Fit2nbbLV& model) {
    
    std::cout << "Summary : Calculate p-value. ";
    // get loglikelihood after marginalization
    auto bestpar = model.GetBestFitParameters();
    double logprob0 = model.LogLikelihood(bestpar);

    // get the fitted function
    auto mean = model.GetFittedFnc();
    // find the 'lower' fitted function
    int nbins = model.GetNbins();
    int downBin = model.GetDownBin();
    int upBin = model.GetUpBin();
    std::vector<int> lCounts(nbins);
    for ( int i = 0; i < nbins; ++i ) lCounts[i] = (int)mean[i];
    
    // counter for p-value
    long int pv = 0;
    long int Niter = 1E05;
    double logprob;

    // random generator
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<> distr(0,1);
    
    ProgressBar bar(Niter);
    // markov chain
    for ( int i = 0; i < Niter; ++i ) {
        bar.Update(i);
        logprob = 0.;
        for ( int j = downBin; j <= upBin; ++j ) {
            
            // propose uniformely +-1 count, then accept with the ratio probability
            if ( distr(eng) < 0.5 ) { if ( mean[j]/(lCounts[j]+1) < distr(eng) ) lCounts[j]++; }

            else                    { if ( mean[j]/(lCounts[j]-1) < distr(eng) ) lCounts[j]--; }
            
            // update loglikelihood
            logprob += BCMath::LogPoisson(lCounts[j], mean[j]);
        } 
        // test
        if ( logprob < logprob0 ) pv++;
    }

    return pv/Niter;
}
