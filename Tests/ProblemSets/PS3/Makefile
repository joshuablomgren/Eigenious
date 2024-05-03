# Makefile for Problem Set 3: Extended Huckel Program

# Compiler
CPP = g++

# Compiler flags
CPPFLAGS = -O2 -std=c++11 

# The directory where the include files are
INC_DIR = ../../../Include

# The directory where the object files are
LIB_DIR = ../../../Lib

# The directory where the executables are
BIN_DIR = ../../../Bin

# All executables
TESTS = ps3_test ps3_energy_diff

factorial: factorial.h
	$(CPP) $(CPPFLAGS) -c factorial.cpp

molecule: molecule.h
	$(CPP) $(CPPFLAGS) -c molecule.cpp 

eht_matrices: eht_matrices.h
	$(CPP) $(CPPFLAGS) -c eht_matrices.cpp

ps3_test: factorial molecule eht_matrices
	$(CPP) $(CPPFLAGS) -o ps3_test test.cpp factorial.o molecule.o eht_matrices.o \
	-I$(INC_DIR) $(LIB_DIR)/eigenious.a 
	mv ps3_test $(BIN_DIR)

ps3_energy_diff: factorial molecule eht_matrices
	$(CPP) $(CPPFLAGS) -o ps3_energy_diff energy_diff.cpp factorial.o molecule.o eht_matrices.o \
	-I$(INC_DIR) $(LIB_DIR)/eigenious.a
	mv ps3_energy_diff $(BIN_DIR)
	rm -f *.o

all: $(TESTS)

clean:
	rm -f *.o 
	cd $(BIN_DIR); rm $(TESTS)

	