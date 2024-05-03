# Eigenious source makefile
# Author: Joshua Blomgren
# Date: December 8, 2023
# Description: This makefile compiles the source code for the Eigenious linear algebra library.
#              It also creates a static library from the compiled object files.


# Compiler
CPP = g++

# Compiler flags
CPPFLAGS = -O2 -std=c++17 -c

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are going to moved afer creation
LIB_DIR     = ../Lib

# Source files
EIGENIOUS_OBJS = matrix.o vector.o operations.o

#Library file
LIBRARY = $(LIB_DIR)/eigenious.a

# Create the Eigenious static library
eigenious.a: $(EIGENIOUS_OBJS)
	ar rcs $(LIBRARY) $(EIGENIOUS_OBJS)
	rm $(EIGENIOUS_OBJS)
	@echo "Eigenious Library created successfully."

# Object files
matrix.o: matrix.cpp $(INC_DIR)/matrix.h 
	$(CPP) $(CPPFLAGS) matrix.cpp -I$(INC_DIR)

vector.o: vector.cpp $(INC_DIR)/vector.h
	$(CPP) $(CPPFLAGS) vector.cpp -I$(INC_DIR)

operations.o: operations.cpp $(INC_DIR)/matrix.h
	$(CPP) $(CPPFLAGS) operations.cpp -I$(INC_DIR)

all: eigenious.a

clean:
	cd $(LIB_DIR); rm -f $(LIBRARY)
