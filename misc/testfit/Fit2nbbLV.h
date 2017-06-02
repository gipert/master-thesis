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

    // GET
    // binning
    std::vector<double> GetBinning() const { return dbin; }
    int GetNbins() { return dbin.size()-1; }
    int GetDownBin() { return downBin; }
    int GetUpBin() { return upBin; }
    // misc
    double Getn2n1() const { return n2n1; }

    // SET
    // initialize ubin vector
    void SetBinning(std::vector<double>& v);
    // set fit boundaries, a=b equals no boundaries
    void SetFitRange(double down, double up);
    // initialize dataBEGe vector
    void SetDataBEGe(std::vector<int>& v) { dataBEGe = v; }
    // initialize simBEGe vector
    void SetSimBEGe(std::vector<std::vector<double>>& v) { simBEGe = v; }

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
    const double n2n1 = 134.594580;
};

#endif
