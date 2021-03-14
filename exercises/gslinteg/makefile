CFLGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY: all
all: out.valueExerciseA.txt

out.valueExerciseA.txt: main Makefile
	./$<

.PHONEY: clean
clean:
	$(RM) main out* *.png