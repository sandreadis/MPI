#
#

# C compiler
CC = gcc
CFLAGS = -g -O3
LIBS = -lm
SRC = serial_mpi.c
TARGET = ./bin/
EXE = serial_mpi

RM_TARGET = rm -f $(TARGET)* 
RM_PRINT = rm -f database/*.txt

all:
		$(CC) $(CFLAGS) $(SRC) -o $(TARGET)$(EXE) $(LIBS)


.PHONY: clean

clean:
		$(RM_TARGET) 
		$(RM_PRINT) 
