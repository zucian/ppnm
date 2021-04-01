CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: psiPlot.png compareplot.png

compareplot.png: compareplot.gpi out.jacobi_timer.txt
	gnuplot $<
psiPlot.png: psiPlot.gpi out.psiPlot.txt
	gnuplot $<

out.jacobi_timer.txt out.psiPlot.txt: main
	./$<

clean:
	$(RM) out.* *.o main *.png