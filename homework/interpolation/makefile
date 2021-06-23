CFLAGS 		 +=  $(shell gsl-config --cflags)
LDLIBS 		 +=  $(shell gsl-config --libs)
LDLIBS 		 +=  -lm -Wall

CC		  		=  gcc
PIPE     		=  cat

INPUT   		=  cosData.txt
LINOUT   		=  linearSplineOutput.txt
QUADOUT  		=  quadraticSplineOutput.txt
CUBICOUT  	=  cubicSplineOutput.txt
EXE   	  	=  main
OBJ  	  		=  main.o linearSpline.o integrateFunction.o quadraticSpline.o cubicSpline.o utilities.o
HDRS    	  =  linearSpline.h integrateFunction.h quadraticSpline.h cubicSpline.h utilities.h
IMGFLS  	  =  linearSplinePlot.png quadraticSplinePlot.png cubicSplinePlot.png
TXT     		=  output.txt

.PHONEY: default
default: $(IMGFLS) print

print: $(TXT)
	$(PIPE) $<

linearSplinePlot.png: linearSplinePlot.gpi linearSplineOutput.txt cosData.txt
	pyxplot $<

quadraticSplinePlot.png: quadraticSplinePlot.gpi quadraticSplineOutput.txt cosData.txt
	pyxplot $<

cubicSplinePlot.png: cubicSplinePlot.gpi cubicSplineOutput.txt cosData.txt
	pyxplot $<

$(LINOUT) $(QUADOUT) $(CUBICOUT) $(TXT): $(EXE) $(INPUT)
	./$< $(INPUT) $(LINOUT) $(QUADOUT) $(CUBICOUT) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(OUTPUT) $(OBJ) $(IMGFLS) $(TXT)