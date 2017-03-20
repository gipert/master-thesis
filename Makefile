# Makefile for GerdaCPT
#
# Author: Luigi Pertoldi (luigi.pertoldi@pd.infn.it)
# Created: 01/02/2017

CC     = c++ -std=c++0x
CFLAGS = $(shell root-config --cflags) \
         $(shell gelatio-config --cflags) \
         $(shell mgdo-config --cflags) \
         $(shell gerda-ada-config --cflags)/gerda-ada -Wall -O3 -g -I./datareader -I./progressbar
LIBS   = $(shell root-config --libs) -lTreePlayer \
         $(shell gelatio-config --libs) \
         $(shell mgdo-config --libs) \
         $(shell gerda-ada-config --libs) -L./lib

all : init lib/libDataReader.so lib/libProgressBar.so bin/getspectra bin/sumspectra

init : 
	-mkdir -p lib bin

lib/libProgressBar.so : progressbar/ProgressBar.cc progressbar/ProgressBar.h
	$(CC) -fPIC -shared -o $@ $<

lib/libDataReader.so : datareader/DataReader.cxx datareader/DataReader.h lib/libProgressBar.so
	$(CC) -fPIC -shared $(CFLAGS) -o $@ $< $(LIBS) -lProgressBar

bin/getspectra : getspectra.cxx lib/libDataReader.so
	$(CC) $(CFLAGS) -o $@ $< $(LIBS) -lDataReader

bin/sumspectra : sumspectra.cxx lib/libDataReader.so
	$(CC) $(CFLAGS) -o $@ $< $(LIBS) -lDataReader -lProgressBar

.PHONY : clean
clean :
	-rm -rf lib bin
