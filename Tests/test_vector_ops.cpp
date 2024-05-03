/**
 * @file test_vector_ops.cpp
 * @brief Test the linear algebra operations on the Vector Class.
 * 
 * This file contains the test cases for performing linear algebra operations on the
 * Vector class.
 * 
 * @author Joshua Blomgren
 * @date December 11, 2023
 */

#include "vector.h"
#include "matrix.h"
#include <cassert>

/**
 * @brief Test the addition of two vectors.
 */
void test_addition() {
    cout << "Testing the addition of two vectors..." << endl;
    Vector<int> vector1 = {1, 2, 3};
    Vector<int> vector2 = {3, 2, 1};
    Vector<int> vector3 = {4, 4, 4};
    assert(vector1 + vector2 == vector3);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the subtraction of two vectors.
 */
void test_subtraction() {
    cout << "Testing the subtraction of two vectors..." << endl;
    Vector<int> vector1 = {4, 4, 4};
    Vector<int> vector2 = {3, 2, 1};
    Vector<int> vector3 = {1, 2, 3};
    assert(vector1 - vector2 == vector3);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the multiplication of a vector and a scalar.
 */
void test_scalar_multiplication() {
    cout << "Testing the multiplication of a vector and a scalar..." << endl;
    Vector<int> vector1 = {1, 2, 3};
    Vector<int> vector2 = {2, 4, 6};
    assert(vector1 * 2 == vector2);
    assert(2 * vector1 == vector2);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test vector multiplication (dot product).
 */
void test_vector_multiplication() {
    cout << "Testing the multiplication of two vectors..." << endl;
    Vector<int> vector1 = {2, 7, 1};
    Vector<int> vector2 = {8, 2, 8};
    assert(vector1 * vector2 == 38);
    assert(vector1.dot(vector2) == 38);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test element-wise multiplication of two vectors.
 */
void test_element_wise_multiplication() {
    cout << "Testing the element-wise multiplication of two vectors..." << endl;
    Vector<int> vector1 = {1, 2, 3};
    Vector<int> vector2 = {3, 2, 1};
    Vector<int> vector3 = {3, 4, 3};
    assert(vector1 % vector2 == vector3);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test division of a vector by a scalar.
 */
void test_scalar_division() {
    cout << "Testing the division of a vector and a scalar..." << endl;
    Vector<int> vector1 = {2, 4, 6};
    Vector<int> vector2 = {1, 2, 3};
    assert(vector1 / 2 == vector2);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test division of two vectors.
 */
void test_vector_division() {
    cout << "Testing the division of two vectors..." << endl;
    Vector<int> vector1 = {2, 4, 6};
    Vector<int> vector2 = {2, 2, 2};
    Vector<int> vector3 = {1, 2, 3};
    assert(vector1 / vector2 == vector3);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the transpose of a vector.
*/
void test_vector_transpose() {
    cout << "Testing the transpose of a vector..." << endl;
    Vector<int> vector1 = {1, 2, 3};
    Matrix<int> matrix1 = {{1}, {2}, {3}};
    assert(vector1.transpose() == matrix1);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test creating a diagonal matrix from a vector.
 */
void test_diagonal_matrix() {
    cout << "Testing creating a diagonal matrix from a vector..." << endl;
    Vector<int> vector1 = {1, 2, 3};
    Matrix<int> matrix1 = {{1, 0, 0}, {0, 2, 0}, {0, 0, 3}};
    assert(diagonal_matrix(vector1) == matrix1);
    cout << "Test passed!" << endl;
}

int main() {
    test_addition();
    test_subtraction();
    test_scalar_multiplication();
    test_vector_multiplication();
    test_element_wise_multiplication();
    test_scalar_division();
    test_vector_division();
    test_vector_transpose();
    test_diagonal_matrix();
    // test_norm();
    return 0;
}