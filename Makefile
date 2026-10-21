PAPER_NAME := Litt_Problem_13_G2

.PHONY: all paper clean

all: paper

paper:
	mkdir -p build output/pdf
	latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build $(PAPER_NAME)_complete.tex
	cp build/$(PAPER_NAME)_complete.pdf output/pdf/$(PAPER_NAME).pdf
	sha256sum output/pdf/$(PAPER_NAME).pdf $(PAPER_NAME)_complete.tex > SHA256SUMS.txt

clean:
	latexmk -C -outdir=build $(PAPER_NAME)_complete.tex
