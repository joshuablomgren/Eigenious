/**
 * @file matrix.cpp
 * @brief Matrix class for linear algebra operations.
 * 
 * The Matrix class allows for the creation of matrices of which matrix operations
 * such as addition, subtraction, multiplication, division can be performed on.
 * The Matrix class is a templated class that allows for the creation of matrices
 * of any type.
 * 
 * @author Joshua Blomgren
 * @date December 8, 2023
 */

#include "matrix.h"

/**
 * @brief Default constructor to create an empty matrix.
 */
template <typename T>
Matrix<T>::Matrix() : n_rows(0), n_cols(0), n_elements(0) {}

/**
 * @brief Constructor to create a matrix of a given size.
 * @param n_rows The number of rows in the matrix.
 * @param n_cols The number of columns in the matrix.
 */
template <typename T>
Matrix<T>::Matrix(int n_rows, int n_cols) : n_rows(n_rows), n_cols(n_cols) {
    if (n_rows < 0 || n_cols < 0) {
        throw invalid_argument("Number of rows and columns must be positive.");
    }
    n_elements = n_rows * n_cols;
    data.resize(n_elements);
}

/**
 * @brief Copy constructor to create a matrix from another matrix.
 * @param matrix The matrix to copy.
 */
template <typename T>
Matrix<T>::Matrix(const Matrix<T>& matrix) {
    // Check that type is int, float, or double.
    if (typeid(T) != typeid(int) && typeid(T) != typeid(float) && typeid(T) != typeid(double)) {
        throw invalid_argument("Matrix type must be int, float, or double.");
    }

    n_rows = matrix.n_rows;
    n_cols = matrix.n_cols;
    n_elements = matrix.n_elements;
    
    // Create new copy of data.
    data.resize(n_elements);
    for (int i = 0; i < n_elements; i++) {
        data[i] = matrix.data[i];
    }
}

/**
 * @brief Constructor to create a matrix from an initializer list of rows.
 * @param list The initializer list of rows.
 * @throw invalid_argument if the matrix type is not int, float, or double.
 * @throw invalid_argument if all rows are not the same size.
 */
template <typename T>
Matrix<T>::Matrix(initializer_list<vector<T>> list) {
    // Check that type is int, float, or double.
    if (typeid(T) != typeid(int) && typeid(T) != typeid(float) && typeid(T) != typeid(double)) {
        throw invalid_argument("Matrix type must be int, float, or double.");
    }

    n_rows = list.size();      
    n_cols = list.begin()-> size(); // Get size of first row using member access operator of iterator.

    // Check that all rows are the same size.
    for (auto row : list) {
        if (row.size() != n_cols) {
            throw invalid_argument("All rows must be the same size.");
        }
    }

    n_elements = n_rows * n_cols;
    data.resize(n_elements);

    // Copy elements from initializer list to matrix.
    int i = 0;
    for (auto row : list) {
        for (auto element : row) {
            data[i] = element;
            i++;
        }
    }
}

/**
 * @brief Accessor to get the element at the given row and column by overloading the
 * () operator.
 * 
 * @param i The row of the element.
 * @param j The column of the element.
 * @return The element at the given row and column.
 */
template <typename T>
T& Matrix<T>::operator()(int row, int col) {
    // Stop program if row or column is out of bounds.
    if (row >= n_rows || col >= n_cols) {
        throw out_of_range("Row or column out of range.");
    }
    if (row < 0 || col < 0) {
        throw out_of_range("Row or column must be non-negative.");
    }
    return data[row * n_cols + col];
}

/**
 * @brief Const accessor to get the element at the given row and column by overloading
 * the () operator.
 * 
 * @param i The row of the element.
 * @param j The column of the element.
 * @return The element at the given row and column.
 */
template <typename T>
const T& Matrix<T>::operator()(int row, int col) const {
    if (row >= n_rows || col >= n_cols) {
        throw out_of_range("Row or column out of range.");
    }
    if (row < 0 || col < 0) {
        throw out_of_range("Row or column must be non-negative.");
    }
    return data[row * n_cols + col];
}

/**
 * @brief Get the number of rows in the matrix.
 * @return The number of rows in the matrix.
 */
template <typename T>
int Matrix<T>::nRows() const {
    return n_rows;
}

/**
 * @brief Get the number of columns in the matrix.
 * @return The number of columns in the matrix.
 */
template <typename T>
int Matrix<T>::nCols() const {
    return n_cols;
}

/**
 * @brief Get the number of elements in the matrix.
 * @return The number of elements in the matrix.
 */
template <typename T>
int Matrix<T>::nElements() const {
    return n_elements;
}

/**
 * @brief Get the size of the matrix (rows, columns).
 * @return The size of the matrix (rows, columns).
 */
template <typename T>
pair <int, int> Matrix<T>::size() const {
    return {n_rows, n_cols};
}

/** 
 * @brief Get a given row of the matrix.
 * @param row The row to get.
 * @return The row of the matrix.
*/
template <typename T>
Vector<T> Matrix<T>::getRow(int row) const {
    if (row < 0 || row >= n_rows) {
        throw out_of_range("Row out of range.");
    }
    Vector<T> row_vector(n_cols);
    for (int i = 0; i < n_cols; i++) {
        row_vector(i) = data[row * n_cols + i];
    }
    return row_vector;
}

/**
 * @brief Get a given column of the matrix.
 * @param col The column to get.
 * @return The column of the matrix.
 */
template <typename T>
Vector<T> Matrix<T>::getCol(int col) const {
    if (col < 0 || col >= n_cols) {
        throw out_of_range("Column out of range.");
    }
    Vector<T> col_vector(n_rows);
    for (int i = 0; i < n_rows; i++) {
        col_vector(i) = data[i * n_cols + col];
    }
    return col_vector;
}

/**
 * @brief Print the matrix.
 */
template <typename T>
void Matrix<T>::print() const {
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            cout << setw(9) << fixed << setprecision(4) << data[i * n_cols + j];
        } 
        cout << endl;
    }
}

/**
 * @brief Resize the matrix while keeping the data.
 * @param n_rows The new number of rows in the matrix.
 * @param n_cols The new number of columns in the matrix.
 */
template <typename T>
void Matrix<T>::resize(int num_rows, int num_cols) {
    // Check if old matrix is bigger than new matrix.
    if (num_rows * num_cols < n_elements) {
        cout << "WARNING: Resizing to a smaller matrix will result in data loss." << endl;
        // Check if the user wants to continue.
        char response;
        cout << "Continue? (y/n): ";
        cin >> response;
        if (response == 'n') {
            return;
        }
    }
    n_rows = num_rows;
    n_cols = num_cols;
    n_elements = n_rows * n_cols;
    data.resize(n_elements);
}

/**
 * @brief Reset matrix with given size.
 * @param n_rows The new number of rows in the matrix.
 * @param n_cols The new number of columns in the matrix.
 */
template <typename T>
void Matrix<T>::reset_size(int num_rows, int num_cols) {
    if (num_rows < 0 || num_cols < 0) {
        throw invalid_argument("Number of rows and columns must be non-negative.");
    }
    n_rows = num_rows;
    n_cols = num_cols;
    n_elements = n_rows * n_cols;
    data.clear();
    data.resize(n_elements);
}

/**
 * @brief Create a matrix of zeros given the size.
 * @param n_rows The number of rows in the matrix.
 * @param n_cols The number of columns in the matrix.
 * @return The matrix of zeros.
 */
template <typename T>
Matrix<T> zeros(int num_rows, int num_cols) {
    if (num_rows < 0 || num_cols < 0) {
        throw invalid_argument("Number of rows and columns must be non-negative.");
    }
    Matrix<T> matrix(num_rows, num_cols);
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            matrix(i, j) = 0;
        }
    }
    return matrix;
}

/**
 * @brief Create a matrix of ones given the size.
 * @param n_rows The number of rows in the matrix.
 * @param n_cols The number of columns in the matrix.
 * @return The matrix of ones.
 */
template <typename T>
Matrix<T> ones(int num_rows, int num_cols) {
    if (num_rows < 0 || num_cols < 0) {
        throw invalid_argument("Number of rows and columns must be non-negative.");
    }
    Matrix<T> matrix(num_rows, num_cols);
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            matrix(i, j) = 1;
        }
    }
    return matrix;
}

/**
 * @brief Fill a matrix with a given value.
 * @param n_rows The number of rows in the matrix.
 * @param n_cols The number of columns in the matrix.
 * @param value The value to fill the matrix with.
 * @return The matrix filled with the given value.
 */
template <typename T>
Matrix<T> fill(int num_rows, int num_cols, T value) {
    if (num_rows < 0 || num_cols < 0) {
        throw invalid_argument("Number of rows and columns must be non-negative.");
    }
    Matrix<T> matrix(num_rows, num_cols);
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            matrix(i, j) = value;
        }
    }
    return matrix;
}

/**
 * @brief Create a matrix of random values between [0,1] given the size.
 * @param n_rows The number of rows in the matrix.
 * @param n_cols The number of columns in the matrix.
 * @return The matrix of random values.
 */
Matrix<double> rand(int num_rows, int num_cols) {
    if (num_rows < 0 || num_cols < 0) {
        throw invalid_argument("Number of rows and columns must be non-negative.");
    }

    // Seed random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(0.0, 1.0);

    Matrix<double> matrix(num_rows, num_cols);
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            matrix(i, j) = dis(gen);
        }
    }
    return matrix;
}

/**
 * @brief Get the absolute value of the matrix.
 * @return The absolute value of the matrix.
 */
template <typename T>
Matrix<T> abs(const Matrix<T>& matrix) {
    Matrix<T> result(matrix.nRows(), matrix.nCols());
    for (int i = 0; i < matrix.nRows(); i++) {
        for (int j = 0; j < matrix.nCols(); j++) {
            result(i, j) = abs(matrix(i, j));
        }
    }
    return result;
}

/**
 * @brief Get the maximum value of the matrix.
 * @return The maximum value of the matrix.
 */
template <typename T>
T Matrix<T>::max() const {
    T max = data[0];
    for (int i = 1; i < n_elements; i++) {
        if (data[i] > max) {
            max = data[i];
        }
    }
    return max;
}


// ------------------------ Template Instantiations ------------------------ //

template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;

template Matrix<int> zeros(int num_rows, int num_cols);
template Matrix<float> zeros(int num_rows, int num_cols);
template Matrix<double> zeros(int num_rows, int num_cols);

template Matrix<int> ones(int num_rows, int num_cols);
template Matrix<float> ones(int num_rows, int num_cols);
template Matrix<double> ones(int num_rows, int num_cols);

template Matrix<int> fill(int num_rows, int num_cols, int value);
template Matrix<float> fill(int num_rows, int num_cols, float value);
template Matrix<double> fill(int num_rows, int num_cols, double value);

template Matrix<int> abs(const Matrix<int>& matrix);
template Matrix<float> abs(const Matrix<float>& matrix);
template Matrix<double> abs(const Matrix<double>& matrix);
