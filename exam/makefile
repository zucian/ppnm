CFLAGS 		 +=  $(shell gsl-config --cflags)
LDLIBS 		 +=  $(shell gsl-config --libs)
LDLIBS 		 +=  -lm -Wall

CC		  		=  gcc
PIPE     		=  cat

INPUT   		=  cosData.txt
CUBICOUT  	=  cubicSplineOutput.txt
SUBOUT 		= subSplineOutput.txt
SUBOUTJUMP		= subSplineOutputJump.txt
JUMPDATA 	= jumpData.txt
EXE   	  	=  main
OBJ  	  		=  main.o integrateFunction.o quadraticSpline.o cubicSpline.o subSpline.o utilities.o
HDRS    	  =  integrateFunction.h quadraticSpline.h cubicSpline.h subSpline.h utilities.h
IMGFLS  	  =  cubicSplinePlot.png subSplinePlot.png subSplinePlotJump.png
TXT     		=  output.txt

.PHONEY: default
default: $(IMGFLS) print

print: $(TXT)
	$(PIPE) $<

cubicSplinePlot.png: cubicSplinePlot.gpi cubicSplineOutput.txt cosData.txt
	pyxplot $<

subSplinePlot.png: subSplinePlot.gpi subSplineOutput.txt cosData.txt
	pyxplot $<

subSplinePlotJump.png: subSplinePlotJump.gpi subSplineOutputJump.txt jumpData.txt
	pyxplot $<

$(LINOUT) $(QUADOUT) $(CUBICOUT) $(SUBOUT) $(TXT): $(EXE) $(INPUT)
	./$< $(INPUT) $(CUBICOUT) $(SUBOUT) $(SUBOUTJUMP) $(JUMPDATA) 1> $(TXT)

$(EXE): $(OBJ) $(HDRS)

clean:
	$(RM) $(EXE) $(OUTPUT) $(OBJ) $(IMGFLS) $(TXT)