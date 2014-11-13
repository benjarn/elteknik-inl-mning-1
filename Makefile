filename=main

run: clean

pdf:
	pdflatex ${filename}.tex

clean: pdf
	rm -f ${filename}.{ps,log,aux,out,dvi,bbl,blg,fls,toc,fdb_latexmk}
