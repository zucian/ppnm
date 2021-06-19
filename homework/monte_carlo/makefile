CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat
PLTUTIL =  pyxplot

ERRPLOT =  errorScaling.png
ERRDATA =  errorScaling.txt
ERRGPI  =  errorScaling.gpi

EXE   	=  main
OBJ  	=  main.o monteCarlo.o
HDRS    =  monteCarlo.h
TXT     =  output.txt

.PHONEY: default
default: print $(ERRPLOT)

$(ERRPLOT): $(ERRGPI) $(ERRDATA)
	pyxplot $<

print: $(TXT)
	$(PIPE) $<

$(TXT) : $(EXE)
	./$< 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ)