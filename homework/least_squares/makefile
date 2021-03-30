CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: fitplot.png out.fitvalues.txt

fitplot.png: plotsetting.gpi out.dataplot.txt out.fitplot.txt
	gnuplot $<

out.fitdat.txt out.fitfunc.txt out.fitvalues.txt: main
	./$<

main: main.c GS_utilities.o

clean:
	$(RM) out.* *.o main *.png