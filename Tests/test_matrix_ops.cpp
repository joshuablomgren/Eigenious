/**
 * @file test_matrix_ops.cpp
 * @brief Test the matrix operations using the Matrix Class.
 * 
 * This file contains the test cases for performing matrix operations using the
 * Matrix class. 
 * 
 * @author Joshua Blomgren
 * @date December 11, 2023
 */

#include "matrix.h"
#include <cassert>
#include <armadillo>

/**
 * @brief Test the addition of two matrices.
 */
void test_addition() {
    cout << "Testing the addition of two matrices..." << endl;
    Matrix<int> matrix1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Matrix<int> matrix2 = {{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};
    Matrix<int> matrix3 = {{10, 10, 10}, {10, 10, 10}, {10, 10, 10}};
    assert(matrix1 + matrix2 == matrix3);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the subtraction of two matrices.
 */
void test_subtraction() {
    cout << "Testing the subtraction of two matrices..." << endl;
    Matrix<int> matrix1 = {{10, 10, 10}, {10, 10, 10}, {10, 10, 10}};
    Matrix<int> matrix2 = {{9, 8, 7}, {6, 5, 4}, {3, 2, 1}};
    Matrix<int> matrix3 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    assert(matrix1 - matrix2 == matrix3);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the multiplication of a matrix and a scalar.
 */
void test_scalar_multiplication() {
    cout << "Testing the multiplication of a matrix and a scalar..." << endl;
    Matrix<int> matrix1 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    Matrix<int> matrix2 = {{2, 4, 6}, {8, 10, 12}, {14, 16, 18}};
    assert(matrix1 * 2 == matrix2);
    assert(2 * matrix1 == matrix2);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the multiplication of two matrices.
*/
void test_matrix_multiplication() {
    cout << "Testing the multiplication of two matrices..." << endl;
    Matrix<int> matrix1 = {{1, 2, 3}, {4, 5, 6}};
    Matrix<int> matrix2 = {{10, 11}, {20, 21}, {30, 31}};
    Matrix<int> matrix3 = {{140, 146}, {320, 335}};
    assert(matrix1 * matrix2 == matrix3);

    cout << "Testing the multiplication of a matrix and a vector..." << endl;
    Matrix<int> matrix4 = {{2, -1}, {1, 1}};
    Vector<int> vector1 = {1, 2};
    Vector<int> vector2 = {0, 3};
    assert(matrix4 * vector1 == vector2);
    assert(vector1 * matrix4 == vector2);

    cout << "Test passed!" << endl;
}

/**
 * @brief Test element-wise multiplication of two matrices.
*/
void test_element_wise_multiplication() {
    cout << "Testing the element-wise multiplication of two matrices..." << endl;
    Matrix<int> matrix1 = {{1, 2}, {3, 4}};
    Matrix<int> matrix2 = {{4, 8}, {0, 5}};
    Matrix<int> matrix3 = {{4, 16}, {0, 20}};
    assert(matrix1 % matrix2 == matrix3);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the division of a matrix and a scalar.
 */
void test_scalar_division() {
    cout << "Testing the division of a matrix and a scalar..." << endl;
    Matrix<int> matrix1 = {{2, 4, 6}, {8, 10, 12}, {14, 16, 18}};
    Matrix<int> matrix2 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    assert(matrix1 / 2 == matrix2);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the division of two matrices.
 */
void test_matrix_division() {
    cout << "Testing the division of two matrices..." << endl;
    Matrix<int> matrix1 = {{10, 20}, {30, 40}};
    Matrix<int> matrix2 = {{2, 4}, {6, 8}};
    Matrix<int> matrix3 = {{5, 5}, {5, 5}};
    assert(matrix1 / matrix2 == matrix3);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test getting inverse of square root of a vector.
*/
void test_vector_inverse_sqrt() {
    cout << "Testing getting inverse of square root of a vector..." << endl;
    Vector<double> vector1 = {1, 4, 9};
    arma::vec vector2 = {1, 4, 9};
    arma::vec arma_vector = 1.0 / arma::sqrt(vector2);
    Vector<double> vector3 = 1.0 / sqrt(vector1);
    cout << "Armadillo: " << endl;
    cout << arma_vector << endl;
    cout << "Eigenious: " << endl;
    cout << vector3 << endl;
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the power of a matrix.
 */
void test_power() {
    cout << "Testing the power of a matrix..." << endl;
    Matrix<int> matrix1 = {{1, 1, 1}, {2, 2, 2}, {3, 3, 3}};
    Matrix<int> matrix2 = {{6, 6, 6}, {12, 12, 12}, {18, 18, 18}};
    assert(matrix1.power(2) == matrix2);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the transpose of a matrix.
 */
void test_transpose() {
    cout << "Testing the transpose of a matrix..." << endl;
    Matrix<int> matrix1 = {{1, 2, 3}, {4, 5, 6}};
    Matrix<int> matrix2 = {{1, 4}, {2, 5}, {3, 6}};
    assert(matrix1.transpose() == matrix2);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the identity matrix.
*/
void test_identity() {
    cout << "Testing the identity matrix..." << endl;
    Matrix<int> matrix1 = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    assert(identity(3) == matrix1);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test eigen_sym function.
 */
void test_eigsym() {
    cout << "Testing the eigen_sym function..." << endl;

    arma::mat A = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    cout << "Armadillo: " << endl;
    cout << eigval << endl;
    cout << eigvec << endl;

    // Eigenious
    Matrix<int> matrix1 = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6}};
    Vector<double> eigval1;
    Matrix<double> eigvec1;
    cout << "Eigenious: " << endl;
    eigen_sym(eigval1, eigvec1, matrix1);
    cout << eigval1 << endl;
    cout << eigvec1 << endl;
    cout << "Test passed!" << endl;
}

/**
 * @brief Test 2-norm function.
 */
void test_norm() {
    cout << "Testing the norm function..." << endl;
    // Armadillo
    arma::mat A = {{1, 2, 3}, {4, 5, 6}};
    cout << "Armadillo: " << endl;
    cout << arma::norm(A) << endl;

    // Eigenious
    Matrix<int> matrix1 = {{1, 2, 3}, {4, 5, 6}};
    cout << "Eigenious: " << endl;
    cout << norm(matrix1) << endl;
    assert (norm(matrix1) == arma::norm(A));
    cout << "Test passed!" << endl;
}

int main() {
    test_addition();
    test_subtraction();
    test_scalar_multiplication();
    test_matrix_multiplication();
    test_element_wise_multiplication();
    test_scalar_division();
    test_matrix_division();
    test_vector_inverse_sqrt();
    test_power();
    test_transpose();
    test_identity();
    test_eigsym();
    test_norm();
    return 0;
}