/**
 * @file test_matrix.cpp
 * @brief Test initialization of the Matrix Class.
 * 
 * This file contains the test cases for creating the Matrix class where
 * normal and edge cases are tested to ensure proper function.
 * 
 * @author Joshua Blomgren
 * @date December 8, 2023
 */

#include "matrix.h"
#include <cassert>

/**
 * @brief Test the default constructor.
 */
void test_default_constructor() {
    cout << "Testing the default constructor..." << endl;
    Matrix<int> matrix;
    assert(matrix.nRows() == 0);
    assert(matrix.nCols() == 0);
    assert(matrix.nElements() == 0);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the constructor with given size.
 */
void test_constructor_with_size() {
    cout << "Testing the constructor with given size..." << endl;
    Matrix<int> matrix(5, 10);
    assert(matrix.nRows() == 5);
    assert(matrix.nCols() == 10);
    assert(matrix.nElements() == 50);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the copy constructor.
 */
void test_copy_constructor() {
    cout << "Testing the copy constructor..." << endl;
    Matrix<int> matrix(5, 10);
    Matrix<int> matrix_copy(matrix);
    assert(matrix_copy == matrix);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the () operator.
 */
void test_assignment_operator() {
    cout << "Testing the () operator..." << endl;
    Matrix<int> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    assert(matrix(0, 0) == 1);
    assert(matrix(0, 1) == 2);
    assert(matrix(0, 2) == 3);
    assert(matrix(1, 0) == 4);
    assert(matrix(1, 1) == 5);
    assert(matrix(1, 2) == 6);
    assert(matrix(2, 0) == 7);
    assert(matrix(2, 1) == 8);
    assert(matrix(2, 2) == 9);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test printing the matrix.
 */
void test_print() {
    cout << "Testing printing the matrix..." << endl;
    Matrix<int> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    matrix.print();
    cout << "Test passed!" << endl;
}

/**
 * @brief Test overloading the << operator to print the matrix.
 */
void test_overload_print() {
    cout << "Testing overloading the << operator to print the matrix..." << endl;
    Matrix<int> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    cout << matrix;
    cout << "Test passed!" << endl;
}

/**
 * @brief Test getting a row of the matrix.
 */
void test_get_row() {
    cout << "Testing getting a row of the matrix..." << endl;
    Matrix<int> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Vector<int> row = matrix.getRow(1);
    assert(row == Vector<int>({4, 5, 6}));
    cout << "Test passed!" << endl;
}

/**
 * @brief Test getting a column of the matrix.
 */
void test_get_col() {
    cout << "Testing getting a column of the matrix..." << endl;
    Matrix<int> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Vector<int> col = matrix.getCol(1);
    assert(col == Vector<int>({2, 5, 8}));
    cout << "Test passed!" << endl;
}

/**
 * @brief Test size of the matrix.
 */
void test_size() {
    cout << "Testing the size function..." << endl;
    Matrix<int> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    assert(matrix.size() == make_pair(3, 3));
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the resize function.
 */
void test_resize() {
    cout << "Testing the resize function..." << endl;
    Matrix<int> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    matrix.resize(2, 2);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the reset_size function.
 */
void test_reset_size() {
    cout << "Testing the reset_size function..." << endl;
    Matrix<int> matrix = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    matrix.reset_size(2, 2);
    assert(matrix.nRows() == 2);
    assert(matrix.nCols() == 2);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the == operator.
 */
void test_equality_operator() {
    cout << "Testing the == operator..." << endl;
    Matrix<int> matrix1 = {{1, 2, 3}, {4, 5, 6}};
    Matrix<int> matrix2 = {{1, 2, 3}, {4, 5, 6}};
    assert(matrix1 == matrix2);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test fill functions (zeros, ones, fill).
 */
void test_fill() {
    cout << "Testing the zeros function..." << endl;
    Matrix<int> matrix = zeros<int>(3, 3);
    assert(matrix == Matrix<int>({{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}));
    
    cout << "Testing the ones function..." << endl;
    matrix = ones<int>(3, 3);
    assert(matrix == Matrix<int>({{1, 1, 1}, {1, 1, 1}, {1, 1, 1}}));

    cout << "Testing the fill function..." << endl;
    Matrix<int> five_matrix = fill(3, 3, 5);
    assert(five_matrix == Matrix<int>({{5, 5, 5}, {5, 5, 5}, {5, 5, 5}}));
    cout << "Test passed!" << endl;
}

/**
 * @brief Test random matrix generator.
 */
void test_rand() {
    cout << "Testing random matrix generation..." << endl;
    Matrix<double> matrix = rand(3, 3);
    cout << matrix;
    cout << "Test passed!" << endl;
}

int main() {
    test_default_constructor();
    test_constructor_with_size();
    test_copy_constructor();
    test_assignment_operator();
    test_print();
    test_overload_print();
    test_get_row();
    test_get_col();
    test_size();
    test_resize();
    test_reset_size();
    test_equality_operator();
    test_fill();
    test_rand();

    return 0;
}

