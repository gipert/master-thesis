# Makefile for GerdaCPT
#
# Author: Luigi Pertoldi (luigi.pertoldi@pd.infn.it)
# Created: 01/02/2017

CC     = c++-4.9 --std=c++11
CFLAGS = $(shell root-config --cflags) \
         $(shell gelatio-config --cflags) \
         $(shell mgdo-config --cflags) \
         $(shell gerda-ada-config --cflags)/gerda-ada -Wall -O3 -g -I./datareader -I./progressbar
LIBS   = $(shell root-config --libs) -lTreePlayer \
         $(shell gelatio-config --libs) \
         $(shell mgdo-config --libs) \
         $(shell gerda-ada-config --libs) -L./datareader -L./progressbar

all : datareader/libDataReader.so progressbar/libProgressBar.so getspectra sumspectra

progressbar/libProgressBar.so : progressbar/ProgressBar.cc progressbar/ProgressBar.h
	$(CC) -fPIC -shared -o $@ $<

datareader/libDataReader.so : datareader/DataReader.cxx datareader/DataReader.h progressbar/libProgressBar.so
	$(CC) -fPIC -shared $(CFLAGS) -o $@ $< $(LIBS) -lProgressBar

getspectra : getspectra.cxx datareader/libDataReader.so
	$(CC) $(CFLAGS) -o $@ $< $(LIBS) -lDataReader

sumspectra : sumspectra.cxx datareader/libDataReader.so
	$(CC) $(CFLAGS) -o $@ $< $(LIBS) -lDataReader -lProgressBar

.PHONY : clean
clean :
	-rm datareader/libDataReader.so
	-rm progressbar/libProgressBar.so
	-rm getspectra
	-rm sumspectra
