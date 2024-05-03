/**
 * @file vector.h
 * @brief Vector class for linear algebra operations.
 * 
 * The Vector class allows for the creation of vectors and the performance of
 * linear algebra operations such as addition, subtraction, multiplication, and division.
 * Vector class is templated to allow for the creation of vectors of any type.
 * 
 * @author Joshua Blomgren
 * @date December 8, 2023
 */

#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <initializer_list>
#include <random>

using namespace std;

// Forward declaration of the Matrix class.
template <typename T>
class Matrix;

/**
 * @class Vector
 * @brief Vector class for linear algebra operations.
 * 
 * The Vector class allows for the creation of vectors and the performance of
 * vector operations such as addition, subtraction, multiplication, and division.
 * The Vector class is a templated class that allows for the creation of vectors
 * of any type.
 * 
 * @tparam T The type of the vector.
 */
template <typename T>
class Vector {
    private:
        vector<T> data;       // The elements of the vector.
        int n_elements;       // The number of elements in the vector.

    public:
        Vector();
        Vector(int n_elements);
        Vector(const Vector<T>& vector);
        Vector(const initializer_list<T>& list);
        T& operator()(int i);
        const T& operator()(int i) const;
        int size() const;
        void print() const;
        void resize(int num_elements);
        void reset_size(int num_elements);

        T dot(const Vector<T>& vector2) const;
        Matrix<T> transpose() const;
        T accumulate() const;
        T max() const;
        

        // ------------------------ Operator Overloads ------------------------ //
        /**
         * @brief Output operator (<<) to print the vector.
         * @param os 
         * @param vector 
         * @return ostream& 
         */
        friend ostream& operator<<(ostream& os, const Vector<T>& vector) {
            for (int i = 0; i < vector.n_elements; i++) {
                os << setw(9) << fixed << setprecision(4) << vector.data[i];
            }
            os << endl;
            return os;
        }

        /**
         * @brief Assignment operator (=) to assign one vector to another.
         * @param vector 
         * @return Vector<T>& 
         */
        Vector<T>& operator=(const Vector<T>& vector) {
            n_elements = vector.n_elements;
            data = vector.data;
            return *this;
        }

        /**
         * @brief Comparison operator (==) to check if two vectors are equal.
         * @param vector2
         * @return bool
         */
        bool operator==(const Vector<T>& vector2) const {
            if (n_elements != vector2.n_elements) {
                return false;
            }
            for (int i = 0; i < n_elements; i++) {
                if (data[i] != vector2.data[i]) {
                    return false;
                }
            }
            return true;
        }

        /**
         * @brief Comparison operator (!=) to check if two vectors are not equal.
         * @param vector2
         * @return bool
         */
        bool operator!=(const Vector<T>& vector2) const {
            return !(*this == vector2);
        }

        /**
         * @brief Addition operator (+) to add two vectors.
         * @param vector2
         * @return Vector<T>
         * @throws invalid_argument if the vectors are not the same size.
         */
        Vector<T> operator+(const Vector<T>& vector2) const {
            if (n_elements != vector2.n_elements) {
                throw invalid_argument("Vectors must have the same number of elements.");
            }
            Vector<T> result(n_elements);
            for (int i = 0; i < n_elements; i++) {
                result.data[i] = data[i] + vector2.data[i];
            }
            return result;
        }

        /**
         * @brief Subtraction operator (-) to subtract two vectors.
         * @param vector2
         * @return Vector<T>
         * @throws invalid_argument if the vectors are not the same size.
         */
        Vector<T> operator-(const Vector<T>& vector2) const {
            if (n_elements != vector2.n_elements) {
                throw invalid_argument("Vectors must have the same number of elements.");
            }
            Vector<T> result(n_elements);
            for (int i = 0; i < n_elements; i++) {
                result.data[i] = data[i] - vector2.data[i];
            }
            return result;
        }

        /**
         * @brief Multiplication operator (*) to multiply vector with a scalar.
         * @param scalar
         * @return Vector<T>
         */
        Vector<T> operator*(T scalar) const {
            Vector<T> result(n_elements);
            for (int i = 0; i < n_elements; i++) {
                result.data[i] = data[i] * scalar;
            }
            return result;
        }

        /**
         * @brief Multiplication operator (*) to multiply a scalar with a vector.
         * @param scalar
         * @param vector2
         * @return Vector<T>
         */
        friend Vector<T> operator*(T scalar, const Vector<T>& vector2) {
            Vector<T> result(vector2.n_elements);
            for (int i = 0; i < vector2.n_elements; i++) {
                result.data[i] = scalar * vector2.data[i];
            }
            return result;
        }

        /**
         * @brief Multiplication operator (*) to multiply two vectors (dot product)
         * @param vector2
         * @return T
         * @throws invalid_argument if the vectors are not the same size.
         */
        T operator*(const Vector<T>& vector2) const {
            if (n_elements != vector2.n_elements) {
                throw invalid_argument("Vectors must have the same number of elements.");
            }
            T result = 0;
            for (int i = 0; i < n_elements; i++) {
                result += data[i] * vector2.data[i];
            }
            return result;
        }

        /**
         * @brief Element-wise multiplication operator (%) to multiply two vectors.
         * @param vector2
         * @return Vector<T>
         * @throws invalid_argument if the vectors are not the same size.
         */
        Vector<T> operator%(const Vector<T>& vector2) const {
            if (n_elements != vector2.n_elements) {
                throw invalid_argument("Vectors must have the same number of elements.");
            }
            Vector<T> result(n_elements);
            for (int i = 0; i < n_elements; i++) {
                result.data[i] = data[i] * vector2.data[i];
            }
            return result;
        }

        /**
         * @brief Division operator (/) to divide a vector by a scalar.
         * @param scalar
         * @return Vector<T>
         */
        Vector<T> operator/(T scalar) const {
            if (scalar == 0) {
                throw invalid_argument("Cannot divide by zero.");
            }
            Vector<T> result(n_elements);
            for (int i = 0; i < n_elements; i++) {
                result.data[i] = data[i] / scalar;
            }
            return result;
        }

        /**
         * @brief Division operator (/) to do element-wise division of two vectors.
         * @param vector2
         * @return Vector<T>
         * @throws invalid_argument if the vectors are not the same size.
         */
        Vector<T> operator/(const Vector<T>& vector2) const {
            if (n_elements != vector2.n_elements) {
                throw invalid_argument("Vectors must have the same number of elements.");
            }
            Vector<T> result(n_elements);
            for (int i = 0; i < n_elements; i++) {
                if (vector2.data[i] == 0 && data[i] != 0) {
                    throw invalid_argument("Cannot divide by zero.");
                }
                result.data[i] = data[i] / vector2.data[i];
            }
            return result;
        }

        /**
         * @brief Division operator (/) to divide a scalar by a vector.
         * @param scalar
         * @param vector2
         * @return Vector<T>
         * @throws invalid_argument if the vector contains a zero.
         */
        friend Vector<T> operator/(T scalar, const Vector<T>& vector2) {
            Vector<T> result(vector2.n_elements);
            for (int i = 0; i < vector2.n_elements; i++) {
                if (vector2.data[i] == 0) {
                    throw invalid_argument("Cannot divide by zero.");
                }
                result.data[i] = scalar / vector2.data[i];
            }
            return result;
        }
};

// ------------------------ Vector Functions ------------------------ //
template <typename T>
Vector<T> zeros(int num_elements);
template <typename T>
Vector<T> ones(int num_elements);
template <typename T>
Vector<T> fill(int num_elements, T value);

Vector<double> rand(int num_elements);

template <typename T>
double norm(const Vector<T>& vector);

template <typename T>
Matrix<T> diagonal_matrix(const Vector<T>& vector);

template <typename T>
Vector<T> sqrt(const Vector<T>& vector);

template <typename T>
Vector<T> abs(const Vector<T>& vector);
