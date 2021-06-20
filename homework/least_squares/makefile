CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: fitplot.png out.fitvaluesexerciseAandB.txt

fitplot.png: plotsetting.gpi out.dataplot.txt out.fitplot.txt
	gnuplot $<

out.dataplot.txt out.fitplot.txt out.fitvaluesexerciseAandB.txt: main
	./$<

main: main.c

clean:
	$(RM) out.* *.o main *.png