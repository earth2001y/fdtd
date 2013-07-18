.SUFFIXES: .tex .dvi .pdf
TEX=platex
TEXOPT=--interaction=nonstopmode

.tex.dvi:
	$(TEX) $(TEXOPT) $< && $(TEX) $(TEXOPT) $< && $(TEX) $(TEXOPT) $<
.dvi.pdf:
	dvipdfmx $<


all: fdtd.pdf

fdtd.pdf: fdtd.dvi
fdtd.dvi: fdtd.tex

clean:
	rm -f *.aux *.log *.dvi *.bib

