CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

exponentialFunction.pdf: exponentialFunction.tex exponentialPlot.pdf
	pdflatex $<
	pdflatex $<

data.txt:main
	./$< > $@

main:main.o

exponentialPlot.pdf: exponentialPlot.gpi data.txt
	pyxplot $<

exponentialFunction.o: main

clean:
	$(RM) exponentialFunction.pdf *.log *.aux out* main *.o fig*.pdf fig*.tex