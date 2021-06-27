CFLAGS +=  $(shell gsl-config --cflags)
LDLIBS +=  $(shell gsl-config --libs)
LDLIBS +=  -lm -Wall

CC   	=  gcc
PIPE    =  cat
PLTUTIL =  pyxplot

EXE   	=  main
OBJ  	=  main.o rootFindingMethods.o GS_utilities.o rungeKutta.o
HDRS    =  rootFindingMethods.h GS_utilities.h rungeKutta.h
TXT     =  output.txt
IMGFL   =  hydrogenPlot.png convergencePlot.png
HGPI    =  hydrogenPlot.gpi
CGPI    =  convergencePlot.gpi
HDATA   =  hydrogenData.txt
CDATA   =  convergenceData.txt

.PHONEY: default
default: print $(IMGFL)

hydrogenPlot.png: $(HGPI) $(HDATA)
	pyxplot $<

convergencePlot.png: $(CGPI) $(CDATA)
	pyxplot $<

print: $(TXT)
#	$(PIPE) $<

$(TXT) : $(EXE)
	./$< $(HDATA) #1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(TXT) $(OBJ) hydrogenPlot.png convergencePlot.png