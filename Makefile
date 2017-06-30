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
COPT = -pdflua -outdir=log -M -MP -MF log/.deps
DIRS = log log/src

all : abstract $(PROJ).pdf

$(PROJ).pdf : $(PROJ).tex $(DIRS) FORCE_MAKE
	@cd img && make all && cd ..
	@echo "Compiling $@..."
	@$(LC) $(COPT) $<
	@cp log/main.pdf ./main.pdf
#if [ ! -h $@ ]; then ln -s log/$@ ./$@; fi

abstract : abstract/abstract.tex FORCE_MAKE
	@cd abstract && make all && cd ..

$(DIRS) :
	@mkdir -p $(DIRS)

clean :
	-rm -rf log abstract/log img/log

.PHONY : FORCE_MAKE clean
-include log/.deps
