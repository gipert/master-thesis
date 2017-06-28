#!/bin/make
#
PROJ = slides
DIRS = log log/src
LC = xelatex -interaction=batchmode
DEPS = $(wildcard img/*.pdf)

all : $(PROJ).pdf

$(PROJ).pdf : $(PROJ).tex $(DIRS) $(DEPS)
	$(LC)  $<
	cp log/slides.pdf ./slides.pdf

$(DIRS) :
	mkdir -p $(DIRS)

clean :
	-rm -rf $(DIRS)

.PHONY : clean
