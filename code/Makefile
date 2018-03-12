#!/bin/make
#
# Author: Luigi Pertoldi (luigi.pertoldi@pd.infn.it)
# Created: 01/02/2017
#
# NOTE: set CC and add system-related variables in misc/vars.mk

CFLAGS    = $(shell root-config --cflags) \
            $(shell gelatio-config --cflags) \
            $(shell mgdo-config --cflags) \
            $(shell gerda-ada-config --cflags)/gerda-ada \
	    $(shell bat-config --cflags) -Wall -O3 -g -I./management -I./progressbar -L./lib
ROOTLIBS  = $(shell root-config --libs) -lTreePlayer
GERDALIBS = $(shell gelatio-config --libs) \
            $(shell mgdo-config --libs) \
            $(shell gerda-ada-config --libs)
DIRS   = bin lib out
include misc/vars.mk

all : $(DIRS) bin/processData bin/processbb bin/sumbb bin/sumbkgext bin/runfit bin/exposure bin/arrayinfo

fit : $(DIRS) bin/runfit

sumtools : $(DIRS) bin/sumbb bin/sumbkgext

$(DIRS)) : 
	-mkdir -p $(DIRS)

lib/libProgressBar.so : progressbar/ProgressBar.cc progressbar/ProgressBar.h
	$(CC) -fPIC -shared -o $@ $<

lib/libDetectorSet.so : management/DetectorSet.cxx management/DetectorSet.h
	$(CC) -fPIC -shared $(CFLAGS) -o $@ $<

lib/libDataReader.so : management/DataReader.cxx management/DataReader.h lib/libProgressBar.so lib/libDetectorSet.so
	$(CC) -fPIC -shared $(CFLAGS) -o $@ $< $(ROOTLIBS) $(GERDALIBS) -lProgressBar -lDetectorSet

lib/libFit2nbbLV.so : fit/Fit2nbbLV.cxx fit/Fit2nbbLV.h
	$(CC) -fPIC -shared $(CFLAGS) -fopenmp -o $@ $< $(shell bat-config --libs)
# ------------------------------------------------------------------------
bin/processData : processData.cxx lib/libDataReader.so
	$(CC) $(CFLAGS) -o $@ $< $(ROOTLIBS) -lDataReader -lDetectorSet

bin/processbb : processbb.cxx lib/libProgressBar.so
	$(CC) $(CFLAGS) -o $@ $< $(ROOTLIBS) -lProgressBar

bin/sumbb : sumbb.cxx lib/libDataReader.so lib/libDetectorSet.so addResolution.cxx
	$(CC) $(CFLAGS) -o $@ $< addResolution.cxx $(ROOTLIBS) -lDataReader -lDetectorSet

bin/sumbkgext : sumbkgext.cxx lib/libDataReader.so lib/libDetectorSet.so addResolution.cxx
	$(CC) $(CFLAGS) -o $@ $< addResolution.cxx $(ROOTLIBS) -lDataReader -lDetectorSet

bin/runfit : fit/runfit.cxx fit/pvalue.cxx lib/libFit2nbbLV.so lib/libProgressBar.so
	$(CC) $(CFLAGS) -fopenmp -o $@ $< fit/pvalue.cxx $(ROOTLIBS) $(shell bat-config --libs) -lFit2nbbLV -lProgressBar

bin/exposure : misc/exposure.cxx lib/libDataReader.so lib/libDetectorSet.so
	$(CC) $(CFLAGS) -o $@ $< $(ROOTLIBS) -lDataReader -lDetectorSet

bin/arrayinfo : misc/arrayinfo.cxx lib/libDetectorSet.so
	$(CC) $(CFLAGS) -o $@ $< $(ROOTLIBS) -lDetectorSet
# ------------------------------------------------------------------------
rundata : 
	bin/processData >/dev/null && misc/sumallbkgext.sh >/dev/null
runbb :
	bin/processbb --2nbb >/dev/null && bin/sumbb --2nbb
runbbLV :
	bin/processbb --2nbbLV >/dev/null && bin/sumbb --2nbbLV

run : rundata runbb runbbLV
	bin/runfit
	telegram-send "make run: task completed"

.PHONY : clean dataclean
clean :
	-rm -rf lib bin

dataclean : 
	-rm data/sumMaGe*
	-rm data/sumData*
