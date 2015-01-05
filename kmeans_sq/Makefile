PROGRAM = main
SOURCES = main.cpp point.cpp dna_strand.cpp
OBJECTS = $(SOURCES:.cpp=.o)
CC = g++

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -Wall -std=c++0x -o $@ $(OBJECTS)

main.o: main.cpp point.h dna_strand.h kmeans_seq.h
	$(CC) -Wall -std=c++0x -c -o main.o main.cpp

point.o: point.cpp point.h
	$(CC) -Wall -std=c++0x -c -o point.o point.cpp

DNA_strand.o: dna_strand.cpp dna_strand.h
	$(CC) -Wall -std=c++0x -c -o dna_strand.o dna_strand.cpp

clean:
	$(RM) -f $(OBJECTS) $(PROGRAM)
