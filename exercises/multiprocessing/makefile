CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += -pthread
CFLAGS += -fopenmp
LDLIBS = -lm -lgomp -lpthread

OBJECTS = main.o randomNumber.o approximatePi.o
HDRS = randomNumber.h approximatePi.h
EXECS = main

default: $(EXECS)
	./$<

$(EXECS): $(OBJECTS)

$(OBJECTS): $(HDRS)

clean:
	$(RM) $(EXECS) $(OBJECTS)