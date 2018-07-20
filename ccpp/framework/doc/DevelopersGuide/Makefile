#
# Makefile for the CCPP (SCM) Users Guide
#

.PHONY: main.pdf all clean

all: main.pdf

main.pdf: main.tex
	latexmk -f -pdf -pdflatex="pdflatex" -use-make main.tex
	open main.pdf

clean:
	latexmk -CA

