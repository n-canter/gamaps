CC=mpicc
BIN=bin/gaMatt

CFLAGS = -pthread\
	-ggdb\
	-Ithird_party\
	-Iinclude\
	-std=gnu99\
	-fopenmp\
	-O3

MATT = $(wildcard src/matt/*.c) 

MPI = $(wildcard src/mpi/*.c)
GA = $(wildcard src/ga/*.c)

obj = $(MATT:.c=.o) $(MPI:.c=.o) $(GA:.c=.o)

$(shell mkdir -p bin)

gamatt: $(obj)
	$(CC) -o $(BIN) $^ $(CFLAGS) -lm


.PHONY: clean

clean:
	rm -f $(obj) bin/gaMatt



