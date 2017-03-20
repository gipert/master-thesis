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
    void SetBinning(std::vector<int>& v) { ubin = v; }
    // initialize dataBEGe vector
    void SetDataBEGe(std::vector<int>& v) { dataBEGe = v; }
    // initialize dataCOAX vector 
    void SetDataCOAX(std::vector<int>& v) { dataCOAX = v; }
    // initialize simBEGe vector
    void SetbbSimBEGe(std::vector<int>& v) { bbSimBEGe = v; }
    // initialize simCOAX vector
    void SetbbSimCOAX(std::vector<int>& v) { bbSimCOAX = v; }

    private:
    
    // new custom data containers (variable binning)
    // WARNING: integer values
    // upper bin bounds
    std::vector<int> ubin;
    // entries for BEGe
    std::vector<int> dataBEGe;
    // entries for enrCOAX
    std::vector<int> dataCOAX;

    // custom model containters (MaGe output)
    // 2nbb
    std::vector<int> bbSimBEGe;
    std::vector<int> bbSimCOAX;
    
    // bkg
    // std::vector<int> bkgSim;
};

#endif
