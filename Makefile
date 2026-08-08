PAPER_NAME := Litt_Problem_13_G2

.PHONY: all paper verify clean

all: paper

paper:
	cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
	mkdir -p output/pdf
	cp paper/main.pdf output/pdf/$(PAPER_NAME).pdf
	sha256sum output/pdf/$(PAPER_NAME).pdf > SHA256SUMS.txt

verify:
	python code/g2_shifted_norm_certificate.py

clean:
	cd paper && latexmk -C main.tex
