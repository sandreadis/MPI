#
#

# C compiler
CC = mpicc
CFLAGS = -g -O3
LIBS = -lm
SRC = knn_mpi.c
TARGET = ./bin/
EXE = knn_mpi

RM_TARGET = rm -f $(TARGET)* 
RM_PRINT = rm -f database/*.txt

all:
		$(CC) $(CFLAGS) $(SRC) -o $(TARGET)$(EXE) $(LIBS)


.PHONY: clean

clean:
		$(RM_TARGET) 
		$(RM_PRINT) 
