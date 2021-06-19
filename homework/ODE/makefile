CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat

OUTPUT  =  harmonicODEOutput.txt SIRODEOutput.txt SIRODEOutput2.txt threeBodyODEOutput.txt
EXE   	=  main
OBJ  	=  main.o rungeKutta.o differentialEquations.o
HDRS    =  rungeKutta.h differentialEquations.h
IMGFL   =  harmonicPlot.png SIRPlot.png SIRPlot2.png threeBodyPlot.png
TXT     =  output.txt

.PHONEY: default
default: print $(IMGFL)

print: $(TXT)
	$(PIPE) $<

harmonicPlot.png: harmonicPlot.gpi $(HARMOUTPUT)
	pyxplot $<

SIRPlot.png: SIRPlot.gpi $(SIROUTPUT)
	pyxplot $<

SIRPlot2.png: SIRPlot2.gpi $(SIR2OUTPUT)
	pyxplot $<

threeBodyPlot.png: threeBodyPlot.gpi $(TPROBOUTPUT)
	pyxplot $<

$(TXT) : $(EXE)
	./$< $(OUTPUT) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OUTPUT) $(OBJ) $(IMGFL)