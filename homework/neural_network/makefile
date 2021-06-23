CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat
PLTUTIL =  pyxplot

EXE   	=  main
OBJ  	=  main.o neuralNetwork.o minimization.o utilities.o inputArray.o
HDRS    =  neuralNetwork.h minimization.h utilities.h inputArray.h
TXT     =  output.txt
DATA    =  cosData.txt networkPrediction.txt
GPI     =  networkPrediction.gpi
IMGFL   =  networkPrediction.png

.PHONEY: default
default: print $(IMGFL)

$(IMGFL): $(GPI) networkPrediction.txt
	pyxplot $<

print: $(TXT)
#	$(PIPE) $<

$(TXT) : $(EXE)
	./$< $(DATA) #1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ) $(IMGFL)