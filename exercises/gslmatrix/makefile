CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm -lgsl -lgslcblas

default: out.txt
	cat out.txt

out.txt: main
	./main 1> out.txt

main: main.o

clean:
	rm main main.o out.txt