CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

all: out.GS_timer.txt

out.GS_timer.txt: main
	./$<

.PHONEY: clean
clean:
	$(RM) out.* *.o main *.png