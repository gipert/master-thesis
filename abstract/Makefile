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
PROJ = abstract
LC   = latexmk
COPT = -silent -pdf -outdir=log -M -MP -MF log/.deps
DIRS = log

all : $(PROJ).pdf

$(PROJ).pdf : $(PROJ).tex $(DIRS) FORCE_MAKE
	@echo "Compiling $@..."
	@$(LC) $(COPT) $<
	@cp log/abstract.pdf ./abstract.pdf
#	if [ ! -h $@ ]; then ln -s log/$@ ./$@; fi

$(DIRS) :
	@mkdir -p $(DIRS)

clean :
	-rm -rf $(DIRS)
	-rm -f $(PROJ).pdf

.PHONY : FORCE_MAKE preview clean
-include log/.deps
