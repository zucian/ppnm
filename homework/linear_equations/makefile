CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

all: out.GS_timer.txt compareplot.png

compareplot.png: compareplot.gpi out.GS_timer.txt
	gnuplot $<

out.GS_timer.txt: main
	./$<

.PHONEY: clean
clean:
	$(RM) out.* *.o main *.png