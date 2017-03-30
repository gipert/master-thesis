/* Fit2nbbLV.h
 *
 * Author: Luigi Pertoldi - luigi.pertoldi@pd.infn.it
 * Created: 20/03/2017
 *
 */

#ifndef __BAT__FIT2NBBLV__H
#define __BAT__FIT2NBBLV__H

#include <BAT/BCModel.h>
#include <vector>
#include <string>

class Fit2nbbLV : public BCModel {
 
    public:
    
    // delete dangerous constructors
    Fit2nbbLV()                            = delete;
    Fit2nbbLV           (Fit2nbbLV const&) = delete;
    Fit2nbbLV& operator=(Fit2nbbLV const&) = delete; 
    // use default destructor
	~Fit2nbbLV()                           = default;

    // custom constructor
	Fit2nbbLV(std::string name = "Fit2nbbLV");

	// methods from BCModel to be overloaded
	double LogLikelihood(const std::vector<double> & parameters);
	//double LogAPrioriProbability(const std::vector<double> & parameters);
    
    // new methods
    // no boundary check for performance reasons
    double Getn2n1() const { return n2n1; }
    double GetBrRatioTl() const { return BrTl; }
    // initialize ubin vector
    void SetBinning(std::vector<double>& v);
    // set fit boundaries, a=b equals no boundaries
    void SetFitRange(double down, double up);
    // initialize dataBEGe vector
    void SetDataBEGe(std::vector<int>& v) { dataBEGe = v; }
    // initialize dataCOAX vector 
    void SetDataCOAX(std::vector<int>& v) { dataCOAX = v; }
    // initialize simBEGe vector
    void SetSimBEGe(std::vector<std::vector<double>>& v) { simBEGe = v; }
    // initialize simCOAX vector
    void SetSimCOAX(std::vector<std::vector<double>>& v) { simCOAX = v; }

    private:
    
    // new custom data containers (variable binning)
    // WARNING: integer values
    // lower bin bounds
    std::vector<double> dbin;
    // upper and lower fit limits
    bool kUseRange;
    int downBin;
    int upBin;
    // entries for BEGe
    std::vector<int> dataBEGe;
    // entries for enrCOAX
    std::vector<int> dataCOAX;

    // custom model containters (MaGe output)
    std::vector<std::vector<double>> simBEGe;
    std::vector<std::vector<double>> simCOAX;
    //
    //// LEGEND:
    //
    // [0] 2bnn
    // [1] 2nbbLV
    // [2] K42homLAr
    // [3] K40onFiberShroud
    // [4] Bi212onFiberShroud
    // [5] Tl208onFiberShroud
    // [6] Pb214onFiberShroud
    // [7] Bi214onFiberShroud
    // [8] alpha
    //

    const double n2n1 = 134.594580;
    const double BrTl = 0.3539;
};

#endif
