# Eigenious test makefile
# Author: Joshua Blomgren
# Date: December 8, 2023
# Description: This makefile creates the executables for the Eigenious tests.

# Compiler
CPP = g++

# Compiler flags
CPPFLAGS = -O2 -std=c++17 -o

# The directory where the include files needed to create the library objects are
INC_DIR = ../Include

# The directory where all object files are moved after creation
LIB_DIR     = ../Lib

# The directory where all test executables are moved after creation
BIN_DIR = ../Bin

# All test executables
TESTS = test_matrix test_vector test_matrix_ops test_vector_ops

# Test files
test_matrix : $(LIB_DIR)/eigenious.a
	$(CPP) $(CPPFLAGS) $(BIN_DIR)/test_matrix test_matrix.cpp \
	-I$(INC_DIR) $(LIB_DIR)/eigenious.a

test_vector : $(LIB_DIR)/eigenious.a
	$(CPP) $(CPPFLAGS) $(BIN_DIR)/test_vector test_vector.cpp \
	-I$(INC_DIR) $(LIB_DIR)/eigenious.a

test_matrix_ops : $(LIB_DIR)/eigenious.a
	$(CPP) $(CPPFLAGS) $(BIN_DIR)/test_matrix_ops test_matrix_ops.cpp \
	-I$(INC_DIR) $(LIB_DIR)/eigenious.a -larmadillo

test_vector_ops : $(LIB_DIR)/eigenious.a
	$(CPP) $(CPPFLAGS) $(BIN_DIR)/test_vector_ops test_vector_ops.cpp \
	-I$(INC_DIR) $(LIB_DIR)/eigenious.a

all: $(TESTS)
	cd ProblemSets/PS3; make all
	cd ProblemSets/PS4; make all

clean:
	cd $(BIN_DIR); rm $(TESTS)
	cd ProblemSets/PS3; make clean
	cd ProblemSets/PS4; make clean