CFLAGS = -o -w -std=gnu99
CFLAGS += $(gsl-config --cflags)
LDLIBS += $(gsl-config --libs)
LDLIBS += -lm -lgsl

all: output ; cat $<
output: main ; ./$< > $@
clean: ; $RM main output