CC = g++
CCFLAGS = -O2
LINALGLIB = -llapack -lgsl -lgslcblas -lz
OMPLIB = -fopenmp

all: DUDE-Seq-2

DUDE-Seq-2: DUDE-Seq-2.cpp
	$(CC) $(CCFLAGS) $< -o $@ $(LINALGLIB) $(OMPLIB)

clean:
	rm -f DUDE-Seq-2
