# main Makefile
#
# Author: Luigi Pertoldi (luigi.pertoldi@pd.infn.it)
# Created: 01/02/2017

CC     = c++-4.9
CFLAGS = $(shell root-config --cflags) -Wall -I./datareader
LIBS   = $(shell root-config --libs) -L./datareader

all : datareader/libDataReader.so main

datareader/libDataReader.so : datareader/DataReader.cxx datareader/DataReader.h
	$(CC) -fPIC -shared $(CFLAGS) -o $@ $< $(LIBS)

main : main.cxx datareader/libDataReader.so
	$(CC) $(CFLAGS) -o $@ $< $(LIBS) -lDataReader

.PHONY : clean
clean :
	-rm datareader/libDataReader.so
	-rm main
