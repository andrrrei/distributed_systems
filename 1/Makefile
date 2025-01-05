CC = mpicc
CFLAGS = -Wall -O2

OBJS = main.o topology.o data.o reduction.o

all: main

main: $(OBJS)
	$(CC) $(CFLAGS) -o main $(OBJS)

main.o: main.c topology.h data.h reduction.h
topology.o: topology.c topology.h
data.o: data.c data.h topology.h
reduction.o: reduction.c reduction.h topology.h

clean:
	rm -f *.o main
