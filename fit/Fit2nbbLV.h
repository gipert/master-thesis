// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__FIT2NBBLV__H
#define __BAT__FIT2NBBLV__H

#include <BAT/BCModel.h>

// This is a Fit2nbbLV header file.
// Model source code is located in file fit/Fit2nbbLV.cxx

// ---------------------------------------------------------
class Fit2nbbLV : public BCModel {
 public:

	// Constructor and destructor
	Fit2nbbLV(const char * name = "Fit2nbbLV");
	~Fit2nbbLV();

	// Methods to overload, see file Fit2nbbLV.cxx
	double LogLikelihood(const std::vector<double> & parameters);
	// double LogAPrioriProbability(const std::vector<double> & parameters);

};
// ---------------------------------------------------------

#endif
