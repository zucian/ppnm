CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

all: out.GS_test.txt

out.GS_test.txt: main
	./$<

.PHONEY: clean
clean:
	$(RM) out.* *.o main *.png