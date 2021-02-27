CC=gcc
CFLAGS+= -Wall -fwrapv -O
LDLIBS=-lm 

.PHONEY: all clean

all: out.txt
	cat $<

out.txt: epsilon makefile
	./epsilon > $@

epsilon: epsilon.o equals.o

epsilon.o: epsilon.c

equal.o: equals.c 

clean:
	$(RM)  epsilon out.txt epsilon.o equals.o