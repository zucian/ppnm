CC		  =  gcc
LDLIBS  =  -lm
ARGS    =  1 2 3 4 5
INPUT   =  input.txt
PIPE    =  cat
TXTS    =  commandRead.txt commandReadError.txt standardInputRead.txt standardInputReadError.txt fileRead.txt fileReadError.txt
EXECS   =  commandRead standardInputRead fileRead

default: $(TXTS)
	$(PIPE) commandRead.txt
	$(PIPE) commandReadError.txt
	$(PIPE) standardInputRead.txt
	$(PIPE) standardInputReadError.txt
	$(PIPE) fileRead.txt
	$(PIPE) fileReadError.txt


$(TXTS): $(EXECS) $(INPUT)
	./commandRead  $$(cat input.txt)  1> commandRead.txt  2> commandReadError.txt
	$(PIPE) $(INPUT) | ./standardInputRead 1> standardInputRead.txt 2> standardInputReadError.txt
	./fileRead $(INPUT) fileRead.txt 			  				2> fileReadError.txt


clean:
	$(RM) $(EXECS) $(TXTS)