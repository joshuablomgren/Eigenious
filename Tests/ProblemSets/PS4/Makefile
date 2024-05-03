# Makefile for Problem Set 4: Build your own SCF program

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
TESTS = ps4_test ps4_oxygen

factorial: factorial.h
	$(CPP) $(CPPFLAGS) -c factorial.cpp

molecule: molecule.h
	$(CPP) $(CPPFLAGS) -c molecule.cpp 

CNDO: CNDO.h
	$(CPP) $(CPPFLAGS) -c CNDO.cpp

ps4_test: factorial molecule CNDO
	$(CPP) $(CPPFLAGS) -o ps4_test test.cpp factorial.o molecule.o CNDO.o \
	-I$(INC_DIR) $(LIB_DIR)/eigenious.a
	mv ps4_test $(BIN_DIR)

ps4_oxygen : factorial molecule CNDO
	$(CPP) $(CPPFLAGS) -o ps4_oxygen oxygen.cpp factorial.o molecule.o CNDO.o \
	-I$(INC_DIR) $(LIB_DIR)/eigenious.a
	mv ps4_oxygen $(BIN_DIR)
	rm -f *.o

all: $(TESTS)

clean:
	rm -f *.o 
	cd $(BIN_DIR); rm $(TESTS)

	