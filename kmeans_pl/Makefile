PROGRAM = main
SOURCES = main.cpp mpi_dna_kmeans.cpp mpi_point_kmeans.cpp
OBJECTS = $(SOURCES:.cpp=.o)
CC = mpic++

all: $(PROGRAM)

$(PROGRAM): main.cpp mpi_kmeans.h
	$(CC) main.cpp mpi_dna_kmeans.cpp mpi_point_kmeans.cpp -o main

clean:
	$(RM) -f $(OBJECTS) $(PROGRAM)
