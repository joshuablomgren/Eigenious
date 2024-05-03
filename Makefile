# Final Project: Eigenious Linear Algebra Library
# Chem 279: Numerical Algorithms applied to Computational Quantum Chemistry
# Created by: Joshua Blomgren
# Date: December 8, 2023

# This makefile creates the library objects for the Eigenious library as well as
# the executables for the Eigenious tests.

# The library objects are created in the Lib directory
# The test executables are created in the Bin directory

SOURCE_DIR = Source
TEST_DIR = Tests

all:
	cd $(SOURCE_DIR); make all
	cd $(TEST_DIR); make all

clean:
	cd $(SOURCE_DIR); make clean
	cd $(TEST_DIR); make clean
	@echo "Cleanup complete."