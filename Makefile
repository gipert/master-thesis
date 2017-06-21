#!/bin/make
#
# makefile for latex projects based on latexmk.
#
# Author: Luigi Pertoldi
# Created: 19/04/2017
#
# Notes: \input does not work with latexmk
#
PROJ = main
LC   = latexmk
COPT = -pdflua -outdir=log -M -MP -MF log/$*.deps
DIRS = log log/src

-include img/Makefile

all : $(PROJ).pdf

$(PROJ).pdf : $(PROJ).tex $(DIRS) FORCE_MAKE
	$(LC) $(COPT) $<
	if [ ! -h $@ ]; then ln -s log/$@ ./$@; fi

preview : $(DIRS)
	$(LC) $(COPT) -pvc $(PROJ).tex

$(DIRS) :
	mkdir -p $(DIRS)

clean :
	$(LC) -CA
	-rm -rf $(DIRS)

.PHONY : FORCE_MAKE preview clean
-include log/*.deps
