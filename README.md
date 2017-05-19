## CONTENTS ##

* `management/`     : utilities to read and manage GERDA data, DataReader and DetectorSet classes
* `progressbar/`    : simple progress bar for c++ loops (git submodule)
* `misc/`
    * `drawspectra.C`: ROOT macro to draw theoretical 2bb spectra
    * `shorten.C`    : ROOT macro to shrink histograms from 8500 to 7500 keV
    * `help`: contains instructions for main program
    * `runconfiguration_mod.db`: modified `runconfiguration.db` from `gerda-metadata` repo
* `fit/`            : class and runfit.cxx to perform the fit
* `data/`           : data needed for the fit
* `Makefile`        : master makefile to compile all the project
* `processData.cxx` : construct energy spectra from data
* `processbb.cxx`   : produce 40 files (in which only one detector acts as a source) containing
                      energy spectra of each detector, starting from MaGe trees (2nbb and 2nbbLV)
* `sumbb.cxx`       : sum the 40 files into one single file (with 40 histograms) properly scaling 
                      each simulation with live time and mass of the detector. The histograms 
                      are then ready to be used in the fit (2nbb and 2nbbLV)
* `sumbkgext.cxx`   : dedicated to external sources of background, scales the histograms for each
                      detector

# NOTES #

* clone with `--recursive` to include ProgressBar submodule
