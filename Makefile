#!/usr/bin/make
#
# Author: Luigi Pertoldi (luigi.pertoldi@pd.infn.it)
# Created: 01/02/2017
#
# NOTE: set CC and add system-related variables in misc/vars.mk

CFLAGS = $(shell root-config --cflags) \
         $(shell gelatio-config --cflags) \
         $(shell mgdo-config --cflags) \
         $(shell gerda-ada-config --cflags)/gerda-ada \
	 $(shell bat-config --cflags) -Wall -O3 -g -I./datareader -I./progressbar
LIBS   = $(shell root-config --libs) -lTreePlayer \
         $(shell gelatio-config --libs) \
         $(shell mgdo-config --libs) \
         $(shell gerda-ada-config --libs) \
	 $(shell bat-config --libs) -L./lib
DIRS   = bin lib out
include misc/vars.mk

all : $(DIRS) bin/getspectra bin/sumspectra bin/runfit

$(DIRS)) : 
	-mkdir -p $(DIRS)

lib/libProgressBar.so : progressbar/ProgressBar.cc progressbar/ProgressBar.h
	$(CC) -fPIC -shared -o $@ $<

lib/libDataReader.so : datareader/DataReader.cxx datareader/DataReader.h lib/libProgressBar.so
	$(CC) -fPIC -shared $(CFLAGS) -o $@ $< $(LIBS) -lProgressBar

lib/libFit2nbbLV.so : fit/Fit2nbbLV.cxx fit/Fit2nbbLV.h
	$(CC) -fPIC -shared $(CFLAGS) -o $@ $< $(LIBS)
# ------------------------------------------------------------------------
bin/getspectra : getspectra.cxx lib/libDataReader.so
	$(CC) $(CFLAGS) -o $@ $< $(LIBS) -lDataReader

bin/sumspectra : sumspectra.cxx lib/libDataReader.so
	$(CC) $(CFLAGS) -o $@ $< $(LIBS) -lDataReader -lProgressBar

bin/runfit : fit/runfit.cxx lib/libFit2nbbLV.so
	$(CC) $(CFLAGS) -o $@ $< $(LIBS) -lFit2nbbLV

.PHONY : clean
clean :
	-rm -rf lib bin
