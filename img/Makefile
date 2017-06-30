#!/bin/make
#
# supports parallel compilation

SRCS  = $(wildcard *.tex)
PDFS  = $(patsubst %.tex, %.pdf, $(SRCS))
LC    = latexmk
FLAGS = -pdf -silent -outdir=log -M -MP -MF log/$*.deps
DIRS  = log

all : $(PDFS)

%.pdf : %.tex  $(DIRS) FORCE_MAKE
	@echo "Compiling $@..."
	@$(LC) $(FLAGS) $<
	@if [ ! -h $@ ]; then ln -s log/$@ ./$@; fi

$(DIRS) :
	mkdir -p $(DIRS)

clean :
	-rm -f $(PDFS)
	-rm -rf log

.PHONY : FORCE_MAKE clean
-include log/*.deps
