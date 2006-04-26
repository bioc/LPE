
default: LPE.pdf

LPE.tex: LPE.Rnw
	echo 'Sweave("$<")' | R --slave

LPE.pdf: LPE.tex
	pdflatex $<
	pdflatex $<
	rm *.log
	rm *.aux
	rm *.out
	rm *.tex
	rm *.toc


%.pdf: %.gif
	convert -density 600x600 %< $@

