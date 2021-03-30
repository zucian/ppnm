CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: psiPlot.png

psiPlot.png: psiPlot.gpi out.psiPlot.txt
	gnuplot $<

out.psiPlot.txt: main
	./$<

clean:
	$(RM) out.* *.o main *.png