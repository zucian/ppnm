default : out.txt
	cat out.txt

out.txt : math                 # math.txt depends on math
	./$< > $@       # this is how you make out.txt from math
              		       # This target basically runs the program and sends it into a .txt file called out.txt

math : math.o       	       		# math depends on math.o
	$(CC) $< -o $@ $(CFLAGS) $(LDFLAGS)    # Now we link the object file to an output program ('-o') called math.

math.o : math.c      	       # math.o depends on math.c
	$(CC) -c $<          # '-c' flag means compile only and do not link.

.PHONEY: clean
clean:                 	       # this is "clean" target it does not depend on anything
	rm -f out.txt math math.o