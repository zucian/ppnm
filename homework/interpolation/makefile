CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: clean all

all: linearSplinePlot.png

linearSplinePlot.png: linearPlotSetting.gpi out.sineData.txt out.linearSpline.txt
	gnuplot $<

out.sineData.txt out.linearSpline.txt: main
	./$<

main: main.c utilities.o linearSpline.o quadraticSpline.o

clean:
	$(RM) out.* *.o main *.png