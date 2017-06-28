#!/bin/make
#
# supports parallel compilation

SRCS  = $(wildcard *.tex)
PDFS  = $(patsubst %.tex, %.pdf, $(SRCS))
LC    = latexmk
FLAGS = -pdf -silent -time -outdir=log -M -MP -MF log/$*.deps
DIRS  = log

all : $(PDFS)

%.pdf : %.tex  $(DIRS) FORCE_MAKE
	$(LC) $(FLAGS) $<
	if [ ! -h $@ ]; then ln -s log/$@ ./$@; fi

$(DIRS) :
	mkdir -p $(DIRS)

clean :
	latexmk $(FLAGS) -C
	rm $(PDFS)
	rm -rf log

.PHONY : FORCE_MAKE clean
-include log/*.deps
