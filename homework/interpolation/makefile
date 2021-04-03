CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

main: main.c utilities.o linearSpline.o

.PHONEY: clean
clean:
	$(RM) out.* *.o main *.png