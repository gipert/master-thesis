#!/bin/make

FILES = $(basename $(wildcard *.tex))
FLAGS = -pdf -pdflatex="pdflatex -interaction=batchmode %O %S" -outdir=log

all :
	latexmk $(FLAGS) -M -MP -MF log/$*.deps;
	for i in $(FILES); do { if [ ! -h $$i.pdf ]; then ln -s log/$$i.pdf ./$$i.pdf; fi; }; done

clean :
	latexmk $(FLAGS) -c

-include log/*.deps
