// test main for GERDA::DataReader class
//
// Author: Luigi Pertoldi (luigi.pertoldi@pd.infn.it)
// Created: 01\02/2017

#include <iostream>

#include "DataReader.h"

int main() {

    GERDA::DataReader reader( "/home/luigi/programs/gerdasw/gerda-metadata",
                              "/home/luigi/code/example/gerda-data",
                              "/home/luigi/code/GERDACPT/runconfiguration_mod.db");

    reader.LoadRun(53);

    
    
    return 0;
}
