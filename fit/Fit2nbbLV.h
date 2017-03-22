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
    // useful in LogLikelihood, no boundary check for performance reasons
    float GetBinCenter(int i) const { return (float)(ubin[i-1]+ubin[i])/2; }
    // initialize ubin vector
    void SetBinning(std::vector<int>& v);
    // set fit boundaries, a=b equals no boundaries
    void SetFitRange(double down, double up);
    // initialize dataBEGe vector
    void SetDataBEGe(std::vector<int>& v) { dataBEGe = v; }
    // initialize dataCOAX vector 
    void SetDataCOAX(std::vector<int>& v) { dataCOAX = v; }
    // initialize simBEGe vector
    void SetSimBEGe(std::vector<std::vector<int>>& v) { simBEGe = v; }
    // initialize simCOAX vector
    void SetSimCOAX(std::vector<std::vector<int>>& v) { simCOAX = v; }

    private:
    
    // new custom data containers (variable binning)
    // WARNING: integer values
    // upper bin bounds
    std::vector<int> ubin;
    // upper and lower fit limits
    bool kUseRange;
    int downBin;
    int upBin;
    // entries for BEGe
    std::vector<int> dataBEGe;
    // entries for enrCOAX
    std::vector<int> dataCOAX;

    // custom model containters (MaGe output)
    // WARNING: must be normalized
    std::vector<std::vector<int>> simBEGe;
    std::vector<std::vector<int>> simCOAX;
    // legend:
    // [0] 2bnn
    // [1] 2nbbLV
    // [2] ...
};

#endif
