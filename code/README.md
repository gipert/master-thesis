# GERDA background model
Under developement.

Contact: luigi.pertoldi@pd.infn.it

## CONTENTS

```text
• management/                 : utilities to read and manage GERDA data, DataReader and DetectorSet classes
• progressbar/                : simple progress bar for c++ loops (git submodule)
• misc/
    • BI.C                    : ROOT macro to calculate background index, works only with 4keV binning
    • exposure.cxx            : program to calculate exposure
    • runconfiguration_mod.db : modified `runconfiguration.db` from `gerda-metadata` repo
    • sumallbkgext.sh         : shell script to reproduce all simulated spectra
    • fixfiles/               : folder with files containing fixing settings for the parameters in the fit
• fit/                        : BAT-derived class and runfit.cxx to perform the fit
• data/                       : data needed for the fit
• processData.cxx             : construct energy spectra from data
• processbb.cxx               : produce 40 files (in which only one detector acts as a source) containing
                                energy spectra of each detector, starting from MaGe trees (2nbb and 2nbbLV)
• sumbb.cxx                   : sum the 40 files into one single file (with 40 histograms) properly scaling 
                                each simulation with live time and mass of the detector. The histograms 
                                are then ready to be used in the fit (2nbb and 2nbbLV)
• sumbkgext.cxx               : dedicated to external sources of background, scales the histograms for each
                                detector
```

### NOTES

* clone with `--recursive` to include ProgressBar submodule
* please check what the Makefile does before doing anything
* set `GERDACPTDIR` env variable to point at the project location!
* just `make fit` if you are only interested in doing the fit with the pre-processed data (required: BAT)
* see `runfit --help`
