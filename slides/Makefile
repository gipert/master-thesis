#!/bin/make
#
# Makefile for latex projects based on latexmk.
# Please set the project name before usage.
# Put source .bib or .tex files in src/
#
# Author: Luigi Pertoldi
# Created: 19/04/2017
#
# Notes: \input does not work with latexmk
#
PROJ = slides
LC   = xelatex
COPT = -shell-escape -interaction nonstopmode -recorder -output-directory=log
DIRS = log log/log

all : $(PROJ).pdf

$(PROJ).pdf : $(PROJ).tex $(DIRS) FORCE_MAKE
	$(LC) $(COPT) $<
	cp log/$(PROJ).pdf $(PROJ).pdf
#	if [ ! -h $@ ]; then ln -s log/$@ ./$@; fi

$(DIRS) :
	mkdir -p $(DIRS)

clean :
	-rm -rf $(DIRS)

.PHONY : clean FORCE_MAKE
