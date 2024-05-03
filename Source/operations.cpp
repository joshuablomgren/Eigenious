/**
 * @file operations.cpp
 * @brief Linear Algebra operations for the Matrix and Vector classes.
 * 
 * This file contains the linear algebra operations able to be performed on the
 * Matrix class and the Vector class that do not require overloading operators.
 * These operations include the power function, transpose, dot product, inverse,
 * and identity.
 * 
 * @author Joshua Blomgren
 * @date December 11, 2023
 */

#include "matrix.h"
#include "vector.h"

/**
 * @brief Vector dot product function.
 * @param vector2 The vector to dot with.
 * @return The dot product of the two vectors.
 */
template <typename T>
T Vector<T>::dot(const Vector<T>& vector2) const {
    if (n_elements != vector2.n_elements) {
        throw invalid_argument("Vectors must be the same size.");
    }
    T result = 0;
    for (int i = 0; i < n_elements; i++) {
        result += data[i] * vector2.data[i];
    }
    return result;
}

/**
 * @brief Matrix power function.
 * @param power The power to raise the matrix to.
 * @return The matrix raised to the given power.
 */
template <typename T>
Matrix<T> Matrix<T>::power(int power) const {
    if (power < 0) {
        throw invalid_argument("Power must be positive.");
    }
    if (n_rows != n_cols) {
        throw invalid_argument("Matrix must be square.");
    }
    Matrix<T> result = *this;
    for (int i = 1; i < power; i++) {
        result = result * *this;
    }
    return result;
}

/**
 * @brief Matrix transpose function.
 * @return The transpose of the matrix.
 */
template <typename T>
Matrix<T> Matrix<T>::transpose() const {
    Matrix<T> result(n_cols, n_rows);
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            result.data[j * n_rows + i] = data[i * n_cols + j];
        }
    }
    return result;
}

/**
 * @brief Vector transpose function.
 * @return The transpose of the vector.
 */
template <typename T>
Matrix<T> Vector<T>::transpose() const {
    Matrix<T> result(n_elements, 1);
    for (int i = 0; i < n_elements; i++) {
        result(i, 0) = data[i];
    }
    return result;
}

/**
 * @brief Create an identity matrix of a given size.
 * @size The size of the identity matrix (square)
*/
Matrix<int> identity(int size) {
    if (size < 0) {
        throw invalid_argument("Size must be non-negative.");
    }
    Matrix<int> result(size, size);
    for (int i = 0; i < size; i++) {
        result(i, i) = 1;
    }
    return result;
}

/**
 * @brief Generate a diagonal matrix from a vector.
 * @param vector The vector to generate the diagonal matrix from.
*/
template <typename T>
Matrix<T> diagonal_matrix(const Vector<T>& vector) {
    Matrix<T> result(vector.size(), vector.size());
    for (int i = 0; i < vector.size(); i++) {
        result(i, i) = vector(i);
    }
    return result;
}

/**
 * @brief Convert T matrix to double matrix.
 * @param matrix The matrix to convert.
 * @return The converted matrix.
*/
template <typename T>
Matrix<double> doubleMatrix(const Matrix<T>& matrix) {
    Matrix<double> result(matrix.nRows(), matrix.nCols());
    for (int i = 0; i < matrix.nRows(); i++) {
        for (int j = 0; j < matrix.nCols(); j++) {
            result(i, j) = static_cast<double>(matrix(i, j));
        }
    }
    return result;
}

/**
 * @brief Check if a matrix is symmetric.
 * @param matrix The matrix to check.
 * @return True if the matrix is symmetric, otherwise false.
*/
template <typename T>
bool isSymmetric(const Matrix<T>& matrix) {
    const double tolerance = 1e-6;

    // Check that the matrix is square.
    if (matrix.nRows() != matrix.nCols()) {
        return false;
    }

    // Check that the matrix is symmetric (within tolerance).
    for (int i = 0; i < matrix.nRows(); i++) {
        for (int j = i + 1; j < matrix.nCols(); j++) {
            if (abs(matrix(i, j) - matrix(j, i)) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Sort inputted eigenvalues and eigenvectors.
 * @param eigenvalues The eigenvalues to sort.
 * @param eigenvectors The eigenvectors to sort.
*/
void eigensort(Vector<double>& eigenvalues, Matrix<double>* eigenvectors = nullptr) {
    int size = eigenvalues.size();
    for (int i = 0; i < size - 1; i++) {
        int min_index = i;
        double min_value = eigenvalues(i);
        for (int j = i; j < size; j++) {
            if (eigenvalues(j) <= min_value) {
                min_index = j;
                min_value = eigenvalues(j);
            }
        }
        if (min_index != i) {
            eigenvalues(min_index) = eigenvalues(i);
            eigenvalues(i) = min_value;

            if (eigenvectors != nullptr) {
                for (int j = 0; j < size; j++) {
                    double temp = (*eigenvectors)(j, i);
                    (*eigenvectors)(j, i) = (*eigenvectors)(j, min_index);
                    (*eigenvectors)(j, min_index) = temp;
                }
            }
        }
    }
}

/**
 * @brief Rotation function needed for Jacobi method.
 * @param matrixA The matrix to rotate.
 * @param s 
 * @param tau
 * @param i
 * @param j
 * @param k
 * @param l
*/
void rotate(Matrix<double>& matrixA, double sin, double tau, int i, int j, int k, int l) {
    double g = matrixA(i, j);
    double h = matrixA(k, l);
    matrixA(i, j) = g - sin * (h + g * tau);
    matrixA(k, l) = h + sin * (g - h * tau);
}

/**
 * @brief Jacobi rotation function.
 * @param matrixA The matrix to rotate.
 * @param i
 * @param j
 * @param eigenvalues
 * @param eigenvectors
 * @param z
 * @param thresh
*/
void jacobiRotate(Matrix<double>& matrixA, int i, int j, Vector<double>& eigenvalues, \
                  Matrix<double>& eigenvectors, Vector<double>& z, double thresh) {
    double offdiag = 100.0 * abs(matrixA(i, j));
    double h = eigenvalues(j) - eigenvalues(i);
    double theta, t, cos, sin, tau;
    double epsilon = numeric_limits<double>::epsilon();

    if (offdiag <= epsilon * abs(eigenvalues(i)) && offdiag <= epsilon * abs(eigenvalues(j))) {
        matrixA(i, j) = 0.0;
    }
    else if (abs(matrixA(i, j)) > thresh) {
        h = eigenvalues(j) - eigenvalues(i);
        if (offdiag <= epsilon * abs(h)) {
            t = matrixA(i, j) / h;
        }
        else {
            theta = 0.5 * h / matrixA(i, j);
            t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0) {
                t = -t;
            }
        }
        cos = 1.0 / sqrt(1.0 + t * t);
        sin = t * cos;
        tau = sin / (1.0 + cos);
        h = t * matrixA(i, j);
        z(i) -= h;
        z(j) += h;
        eigenvalues(i) -= h;
        eigenvalues(j) += h;
        matrixA(i, j) = 0.0;
        for (int k = 0; k < i; k++) {
            rotate(matrixA, sin, tau, k, i, k, j);
        }
        for (int k = i + 1; k < j; k++) {
            rotate(matrixA, sin, tau, i, k, k, j);
        }
        for (int k = j + 1; k < matrixA.nRows(); k++) {
            rotate(matrixA, sin, tau, i, k, j, k);
        }
        for (int k = 0; k < matrixA.nRows(); k++) {
            rotate(eigenvectors, sin, tau, k, i, k, j);
        }
    }
}



/**
 * @brief Solve for the eigenvalues and eigenvectors of a symmetric matrix.
 * 
 * This functions calculates the eigenvalues and eigenvectors of a symmetric
 * matrix using the Jacobi method. The eigenvalues are stored in a vector and
 * the eigenvectors are stored in a matrix.
 * 
 * @param eigenvalues The vector to store the eigenvalues in.
 * @param eigenvectors The matrix to store the eigenvectors in.
 * @param matrixA The matrix to solve for the eigenvalues and eigenvectors of.
 * @throws invalid_argument if the matrix is not symmetric.
 * @throws runtime_error if the matrix does not converge.
 */
template <typename T>
void eigen_sym(Vector<double>& eigenvalues, Matrix<double>& eigenvectors, const Matrix<T>& matrix) {
    // Check if the matrix is symmetric.
    if (!isSymmetric(matrix)) {
        throw invalid_argument("Matrix must be symmetric.");
    }

    // Initialize values and matrices.
    Matrix<double> matrixA = doubleMatrix(matrix);
    const int size = matrix.nRows();
    eigenvalues.reset_size(size);
    eigenvectors.reset_size(size, size);
    int rotations = 0;
    const double epsilon = numeric_limits<double>::epsilon();
    double threshold;

    // Initalize eigenvectors to the identity matrix.
    for (int i = 0; i < size; i++) {
        eigenvectors(i, i) = 1.0;
    }

    // Initialize eigenvalues to the diagonal of the matrix.
    Vector<double> b = Vector<double>(size);
    Vector<double> z = Vector<double>(size);
    for (int i = 0; i < size; i++) {
        eigenvalues(i) = matrixA(i, i);
        b(i) = eigenvalues(i);
        z(i) = 0.0;
    }


    // Jacobi Iteration
    for (int iter = 1; iter <= 50; iter++) {
        double sum = 0.0;

        // Sum off-diagonal elements.
        for (int i = 0; i < size - 1; i++) {
            for (int j = i + 1; j < size; j++) {
                sum += abs(matrixA(i, j));
            }
        }

        // Check for convergence.
        if (sum == 0.0) {
            eigensort(eigenvalues, &eigenvectors);
            return;
        }

        // First three sweeps.
        if (iter < 4) {
            threshold = 0.2 * sum / (size * size);
        }
        else {
            threshold = 0.0;
        }
        
        // Iterate through off-diagonal elements.
        for (int i = 0; i < size - 1; i++) {
            for (int j = i + 1; j < size; j++) {
                double offdiag = 100.0 * abs(matrixA(i, j));

                // After four times, skip the rotation if off-diagonal is small
                if (iter > 4 && offdiag <= epsilon * abs(eigenvalues(i)) && \
                    offdiag <= epsilon * abs(eigenvalues(j))) {
                    matrixA(i, j) = 0.0;
                }
                else if (abs(matrixA(i, j)) > threshold) {
                    // Jacobi Rotation
                    jacobiRotate(matrixA, i, j, eigenvalues, eigenvectors, z, threshold);
                    rotations++;
                }
            }
        }

        // Update eigenvalues and reset z.
        for (int i = 0; i < size; i++) {
            b(i) += z(i);
            eigenvalues(i) = b(i);
            z(i) = 0.0;
        }
    }

    // Throw error if the matrix did not converge.
    throw runtime_error("Too many iterations, Matrix did not converge.");
}

/**
 * @brief Calculate the 2-norm of a vector.
 * @param vector The vector to calculate the norm of.
 * @return The 2-norm of the vector. 
*/
template <typename T>
double norm(const Vector<T>& vector) {
    double result = 0.0;
    for (int i = 0; i < vector.size(); i++) {
        result += vector(i) * vector(i);
    }
    return sqrt(result);
}

/**
 * @brief Calculate the 2-norm of a matrix.
 * @param matrix The matrix to calculate the norm of.
 * @return The 2-norm of the matrix.
*/
template <typename T>
double norm(const Matrix<T>& matrix) {
    // Get eigenvalues of matrix.transpose() * matrix
    Vector<double> eigenvalues;
    Matrix<double> eigenvectors;
    eigen_sym(eigenvalues, eigenvectors, matrix.transpose() * matrix);

    // Find the largest eigenvalue
    double max = 0.0;
    for (int i = 0; i < eigenvalues.size(); i++) {
        if (eigenvalues(i) > max) {
            max = eigenvalues(i);
        }
    }

    // return the square root of the largest eigenvalue
    return sqrt(max);
}


// ------------------------ Functions ------------------------ //
template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;
template class Vector<int>;
template class Vector<float>;
template class Vector<double>;

template Matrix<int> diagonal_matrix(const Vector<int>& vector);
template Matrix<float> diagonal_matrix(const Vector<float>& vector);
template Matrix<double> diagonal_matrix(const Vector<double>& vector);

template Matrix<double> doubleMatrix(const Matrix<int>& matrix);
template Matrix<double> doubleMatrix(const Matrix<float>& matrix);
template Matrix<double> doubleMatrix(const Matrix<double>& matrix);

template void eigen_sym(Vector<double>& eigenvalues, Matrix<double>& eigenvectors, const Matrix<int>& matrix);
template void eigen_sym(Vector<double>& eigenvalues, Matrix<double>& eigenvectors, const Matrix<float>& matrix);
template void eigen_sym(Vector<double>& eigenvalues, Matrix<double>& eigenvectors, const Matrix<double>& matrix);

template double norm(const Vector<int>& vector);
template double norm(const Vector<float>& vector);
template double norm(const Vector<double>& vector);

template double norm(const Matrix<int>& matrix);
template double norm(const Matrix<float>& matrix);
template double norm(const Matrix<double>& matrix);