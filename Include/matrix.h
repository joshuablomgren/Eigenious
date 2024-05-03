/**
 * @file matrix.h   
 * @brief Matrix Class for linear algebra operations.
 * 
 * The Matrix class allows for the creation of matrices of which matrix operations
 * such as addition, subtraction, multiplication, division can be performed on.
 * The Matrix class is a templated class that allows for the creation of matrices
 * of any type.
 * 
 * @author Joshua Blomgren
 * @date December 8, 2023
 */
#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <initializer_list>
#include <random>
#include <limits>
#include "vector.h"

using namespace std;

/**
 * @class Matrix
 * @brief Matrix class for linear algebra operations.
 * 
 * The Matrix class allows for the creation of matrices and the performance of
 * matrix operations such as addition, subtraction, multiplication, and division.
 * The Matrix class is a templated class that allows for the creation of matrices
 * of any type.
 * 
 * @tparam T The type of the matrix.
 */
template <typename T>
class Matrix {
    private:
        vector<T> data;       // The elements of the matrix.
        int n_rows;           // The number of rows in the matrix.
        int n_cols;           // The number of columns in the matrix.
        int n_elements;       // The number of elements in the matrix.

    public:
        Matrix();
        Matrix(int n_rows, int n_cols);
        Matrix(const Matrix<T>& matrix);
        Matrix(initializer_list<vector<T>> list);
        T& operator()(int i, int j);
        const T& operator()(int i, int j) const;
        int nRows() const;
        int nCols() const;
        int nElements() const;
        pair <int, int> size() const;
        Vector<T> getRow(int row) const;
        Vector<T> getCol(int col) const;
        void print() const;
        void resize(int num_rows, int num_cols);
        void reset_size(int num_rows, int num_cols);

        Matrix<T> power(int power) const;
        Matrix<T> transpose() const;
        T max() const;


        // ------------------------ Operator Overloads ------------------------ //
        
        friend ostream& operator<<(ostream& os, const Matrix<T>& matrix) {
            for (int i = 0; i < matrix.n_rows; i++) {
                for (int j = 0; j < matrix.n_cols; j++) {
                    if (matrix.data[i * matrix.n_cols + j] == 0) {
                        os << setw(10) << 0;
                    }
                    else {
                        os << setw(10) << fixed << setprecision(4) << matrix.data[i * matrix.n_cols + j];
                    }
                }
                os << '\n';
            }
            return os;
        }

        /**
         * @brief Assignment operator (=) to assign one matrix to another.
         * @param matrix2
         * @return Matrix<T>& 
         */
        Matrix<T>& operator=(const Matrix<T>& matrix2) {
            n_rows = matrix2.n_rows;
            n_cols = matrix2.n_cols;
            n_elements = matrix2.n_elements;
            data = matrix2.data;
            return *this;
        }

        /**
         * @brief Comparison operator (==) to check if two matrices are equal.
         * @param matrix2 
         * @return true 
         * @return false 
         */
        bool operator==(const Matrix<T>& matrix2) const {
            if (n_rows != matrix2.n_rows || n_cols != matrix2.n_cols) {
                return false;
            }
            for (int i = 0; i < n_elements; i++) {
                if (data[i] != matrix2.data[i]) {
                    return false;
                }
            }
            return true;
        }

        /**
         * @brief Comparison operator (!=) to check if two matrices are not equal.
         * @param matrix2
         * @return true
         * @return false
         */
        bool operator!=(const Matrix<T>& matrix2) const {
            if (n_rows != matrix2.n_rows || n_cols != matrix2.n_cols) {
                return true;
            }
            for (int i = 0; i < n_elements; i++) {
                if (data[i] != matrix2.data[i]) {
                    return true;
                }
            }
            return false;
        }

        /**
         * @brief Addition operator (+) to add two matrices.
         * @param matrix2 
         * @return Matrix<T> 
         */
        Matrix<T> operator+(const Matrix<T>& matrix2) const {
            if (n_rows != matrix2.n_rows || n_cols != matrix2.n_cols) {
                throw invalid_argument("Matrices must be the same size to add.");
            }
            Matrix<T> matrix(n_rows, n_cols);
            for (int i = 0; i < n_elements; i++) {
                matrix.data[i] = data[i] + matrix2.data[i];
            }
            return matrix;
        }

        /**
         * @brief Subtraction operator (-) to subtract two matrices.
         * @param matrix2 
         * @return Matrix<T> 
         */
        Matrix<T> operator-(const Matrix<T>& matrix2) const {
            if (n_rows != matrix2.n_rows || n_cols != matrix2.n_cols) {
                throw invalid_argument("Matrices must be the same size to subtract.");
            }
            Matrix<T> matrix(n_rows, n_cols);
            for (int i = 0; i < n_elements; i++) {
                matrix.data[i] = data[i] - matrix2.data[i];
            }
            return matrix;
        }

        /**
         * @brief Multiplication operator (*) to multiply matrix with a scalar.
         * @param scalar
         * @return Matrix<T> 
         */
        Matrix<T> operator*(const T& scalar) const {
            Matrix<T> matrix(n_rows, n_cols);
            for (int i = 0; i < n_elements; i++) {
                matrix.data[i] = data[i] * scalar;
            }
            return matrix;
        }

        /**
         * @brief Multiplication operator (*) to a scalar with a matrix.
         * @param scalar
         * @param matrix2
         * @return Matrix<T> 
         */
        friend Matrix<T> operator*(const T& scalar, const Matrix<T>& matrix2) {
            Matrix<T> matrix(matrix2.n_rows, matrix2.n_cols);
            for (int i = 0; i < matrix2.n_elements; i++) {
                matrix.data[i] = scalar * matrix2.data[i];
            }
            return matrix;
        }

        /**
         * @brief Multiplication operator (*) to multiply two matrices.
         * @param matrix2
         * @return Matrix<T> 
         */
        Matrix<T> operator*(const Matrix<T>& matrix2) const {
            if (n_cols != matrix2.n_rows) {
                throw invalid_argument("Number of columns in matrix 1 must equal number of rows in matrix 2 to do matrix * matrix.");
            }
            Matrix<T> matrix(n_rows, matrix2.n_cols);
            for (int i = 0; i < n_rows; i++) {
                for (int j = 0; j < matrix2.n_cols; j++) {
                    for (int k = 0; k < n_cols; k++) {
                        matrix(i, j) += data[i * n_cols + k] * matrix2.data[k * matrix2.n_cols + j];
                    }
                }
            }
            return matrix;
        }

        /**
         * @brief Multiplication operator (*) to multiply a matrix with a vector.
         * @param vector2
         * @return Vector<T>
         * @throws invalid_argument if the number of columns in the matrix does not equal vector size. 
         */
        Vector<T> operator*(const Vector<T>& vector2) const {
            if (n_cols != vector2.size()) {
                throw invalid_argument("Number of columns in matrix must be equal to vector size to do matrix * vector.");
            }
            Vector<T> vector(n_rows);
            for (int i = 0; i < n_rows; i++) {
                for (int j = 0; j < n_cols; j++) {
                    vector(i) += data[i * n_cols + j] * vector2(j);
                }
            }
            return vector;
        }

        /**
         * @brief Multiplication operator (*) to multiply a vector with a matrix.
         * @param vector2
         * @param matrix2
         * @return Vector<T>
         * @throws invalid_argument if the number of rows in the matrix does not equal vector size. 
         */
        friend Vector<T> operator*(const Vector<T>& vector2, const Matrix<T>& matrix2) {
            if (matrix2.n_cols != vector2.size()) {
                throw invalid_argument("Number of columns in matrix must be equal to vector size to do vector * matrix.");
            }
            Vector<T> vector(matrix2.n_rows);
            for (int i = 0; i < matrix2.n_rows; i++) {
                for (int j = 0; j < matrix2.n_cols; j++) {
                    vector(i) += vector2(j) * matrix2.data[i * matrix2.n_cols + j];
                }
            }
            return vector;
        }

        /**
         * @brief Multiplication operator (%) to do element-wise multiplication of two matrices.
         * @param matrix2
         * @return Matrix<T> 
         */
        Matrix<T> operator%(const Matrix<T>& matrix2) const {
            if (n_rows != matrix2.n_rows || n_cols != matrix2.n_cols) {
                throw invalid_argument("Matrices must be the same size to do element-wise multiplication.");
            }
            Matrix<T> matrix(n_rows, n_cols);
            for (int i = 0; i < n_elements; i++) {
                matrix.data[i] = data[i] * matrix2.data[i];
            }
            return matrix;
        }

        /**
         * @brief Division operator (/) to divide a matrix by a scalar.
         * @param scalar
         * @return Matrix<T> 
         */
        Matrix<T> operator/(const T& scalar) const {
            if (scalar == 0) {
                throw invalid_argument("Cannot divide by zero.");
            }
            Matrix<T> matrix(n_rows, n_cols);
            for (int i = 0; i < n_elements; i++) {
                matrix.data[i] = data[i] / scalar;
            }
            return matrix;
        }

        /**
         * @brief Division operator (/) to do element-wise division of two matrices.
         * @param matrix2
         * @return Matrix<T> 
         */
        Matrix<T> operator/(const Matrix<T>& matrix2) const {
            if (n_rows != matrix2.n_rows || n_cols != matrix2.n_cols) {
                throw invalid_argument("Matrices must be the same size to do element-wise division.");
            }
            Matrix<T> matrix(n_rows, n_cols);
            for (int i = 0; i < n_elements; i++) {
                if (matrix2.data[i] == 0 && data[i] != 0) {
                    throw invalid_argument("Cannot divide by zero.");
                }
                matrix.data[i] = data[i] / matrix2.data[i];
            }
            return matrix;
        }
};

// ------------------------ Matrix Functions ------------------------ //
template <typename T>
Matrix<T> zeros(int num_rows, int num_cols);
template <typename T>
Matrix<T> ones(int num_rows, int num_cols);
template <typename T>
Matrix<T> fill(int num_rows, int num_cols, T value);
Matrix<double> rand(int num_rows, int num_cols);

Matrix<int> identity(int size);

template <typename T>
Matrix<double> doubleMatrix(const Matrix<T>& matrix);

template <typename T>
bool isSymmetric(const Matrix<T>& matrix);

void rotate(Matrix<double>& matrixA, double s, double tau, int i, int j, int k, int l);
void jacobiRotate(Matrix<double>& matrixA, int i, int j, Vector<double>& eigenvalues, \
                  Matrix<double>& eigenvectors, Vector<double>& z, double thresh);

template <typename T>
void eigen_sym(Vector<double>& eigenvalues, Matrix<double>& eigenvectors, const Matrix<T>& matrix);

template <typename T>
double norm(const Matrix<T>& matrix);

template <typename T>
Matrix<T> abs(const Matrix<T>& matrix);
