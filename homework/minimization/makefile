CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat
PLTUTIL =  pyxplot

EXE   	=  main
OBJ  	=  main.o minimization.o inputArray.o
HDRS    =  minimization.h inputArray.h
TXT     =  output.txt
IMGFL   =  higgsFit.png
GPI     =  higgsFit.gpi
DATA    =  higgsFit.txt

.PHONEY: default
default: print $(IMGFL)

higgsFit.png: $(GPI) $(DATA)
	pyxplot $<

print: $(TXT)
#	$(PIPE) $<

$(TXT) : $(EXE)
	./$< $(DATA) #1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ) higgsFit.png