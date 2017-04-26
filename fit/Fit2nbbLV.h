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
	double LogLikelihood(const std::vector<double>& parameters);
	//double LogAPrioriProbability(const std::vector<double> & parameters);
    
    // GET
    // binning
    std::vector<double> GetBinning() const { return dbin; }
    int GetNbins() { return dbin.size()-1; }
    int GetDownBin() { return downBin; }
    int GetUpBin() { return upBin; }
    // functions
    std::vector<double> GetFittedFncBEGe(std::vector<double>& bestpar);
    std::vector<double> GetFittedFncCOAX(std::vector<double>& bestpar);
    std::vector<int> GetDataBEGe() const { return dataBEGe; }
    std::vector<int> GetDataCOAX() const { return dataCOAX; }
    // misc
    double Getn2n1() const { return n2n1; }
    double GetBrRatioTl() const { return BrTl; }
   
    // SET
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

    // write histos on file
    void WriteHistosOnFile(std::string filename);

    private:
    
    // new custom data containers (variable binning)
    // WARNING: integer values
    // lower bin bounds, size = nbins+1
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
    // [3] K40fibers <-----
    // [4] Bi212fibers
    // [5] Tl208fibers
    // [6] Pb214fibers
    // [7] Bi214fibers
    // [8] alpha
    // [9] nPlus
    // [10] pPlus
    // [11] Ac228holder <----- 
    // [12] Co60holder
    // [13] K40holder
    // [14] Bi212holder 
    // [15] Tl208holder 
    // [16] Pb214holder 
    // [17] Bi214holder 
    // [18] K40cables  <-----
    // [19] Bi212cables
    // [20] Tl208cables
    // [21] Pb214cables
    // [22] Bi214cables
    // [23] K40minishroud <-----
    // [24] Pa234minishroud
    // [25] Bi207minishroud
    // [26] Bi207cables
    // [27] Bi207holder
    // [28] Pb214minishroud
    // [29] Bi214minishroud
    // [30] K42minishroudsurface
    // [31] Pa234cables
    // [32] Pa234holder
    // [33] K42homLArAA

    const double n2n1 = 134.594580;
    const double BrTl = 0.3539;
};

#endif
