/**
 * @file vector.cpp
 * @brief Vector class for linear algebra operations.
 * 
 * The Vector class allows for the creation of vectors and the performance of
 * linear algebra operations such as addition, subtraction, multiplication, and division.
 * Vector class is templated to allow for the creation of vectors of any type.
 * 
 * @author Joshua Blomgren
 * @date December 8, 2023
 */

#include "vector.h"

/**
 * @brief Default constructor to create an empty vector.
 */
template <typename T>
Vector<T>::Vector() : n_elements(0) {}

/**
 * @brief Constructor to create a vector of a given size.
 * @param n_elements The number of elements in the vector.
 */
template <typename T>
Vector<T>::Vector(int n_elements) : n_elements(n_elements) {
    if (n_elements < 0) {
        throw invalid_argument("Number of elements must be non-negative.");
    }
    data.resize(n_elements);
}

/**
 * @brief Copy constructor to create a vector from another vector.
 * @param vector The vector to copy.
 */
template <typename T>
Vector<T>::Vector(const Vector<T>& vector) {
    n_elements = vector.n_elements;
    data = vector.data;
}

/**
 * @brief Constructor to create a vector from an initializer list.
 * @param list The initializer list.
 */
template <typename T>
Vector<T>::Vector(const initializer_list<T>& list) {
    data = vector<T>(list);
    n_elements = data.size();
}

/**
 * @brief Overload the () operator to access elements of the vector.
 * @param i The index of the element to access.
 * @return The element at the given index.
 */
template <typename T>
T& Vector<T>::operator()(int i) {
    if (i < 0 || i >= n_elements) {
        throw invalid_argument("Index out of bounds.");
    }
    return data[i];
}

/**
 * @brief Const overload the () operator to access elements of the vector.
 * @param i The index of the element to access.
 * @return The element at the given index.
 */
template <typename T>
const T& Vector<T>::operator()(int i) const {
    if (i < 0 || i >= n_elements) {
        throw invalid_argument("Index out of bounds.");
    }
    return data[i];
}

/**
 * @brief Get the number of elements in the vector.
 * @return The number of elements in the vector.
 */
template <typename T>
int Vector<T>::size() const {
    return n_elements;
}

/**
 * @brief Print the vector.
 */
template <typename T>
void Vector<T>::print() const {
    for (int i = 0; i < n_elements; i++) {
        cout << setw(9) << fixed << setprecision(4) << data[i];
    }
    cout << endl;
}

/**
 * @brief Resize the vector keeping the current elements.
 * @param num_elements The new number of elements in the vector.
 */
template <typename T>
void Vector<T>::resize(int num_elements) {
    if (num_elements < 0) {
        throw invalid_argument("Number of elements must be non-negative.");
    }
    if (num_elements < n_elements) {
        cout << "WARNING: resizing to a smaller size may result in loss of data." << endl;
        // Check user response
        char response;
        cout << "Continue? (y/n): ";
        cin >> response;
        if (response == 'n') {
            return;
        }
    }
    n_elements = num_elements;
    data.resize(n_elements);
}

/**
 * @brief Reset the matrix with a new size.
 * @param num_elements The new number of elements in the vector.
 */
template <typename T>
void Vector<T>::reset_size(int num_elements) {
    if (num_elements < 0) {
        throw invalid_argument("Number of elements must be non-negative.");
    }
    n_elements = num_elements;
    data.resize(n_elements);
}

/**
 * @brief Create a vector of zeros with a given size.
 * @param num_elements The number of elements in the vector.
 * @return The vector of zeros.
 */
template <typename T>
Vector<T> zeros(int num_elements) {
    if (num_elements < 0) {
        throw invalid_argument("Number of elements must be non-negative.");
    }
    Vector<T> vector(num_elements);
    for (int i = 0; i < num_elements; i++) {
        vector(i) = 0;
    }
    return vector;
}

/**
 * @brief Create a vector of ones with a given size.
 * @param num_elements The number of elements in the vector.
 * @return The vector of ones.
 */
template <typename T>
Vector<T> ones(int num_elements) {
    if (num_elements < 0) {
        throw invalid_argument("Number of elements must be non-negative.");
    }
    Vector<T> vector(num_elements);
    for (int i = 0; i < num_elements; i++) {
        vector(i) = 1;
    }
    return vector;
}

/**
 * @brief Create a vector of a given size with a given value.
 * @param num_elements The number of elements in the vector.
 * @param value The value to fill the vector with.
 * @return The vector of the given value.
 */
template <typename T>
Vector<T> fill(int num_elements, T value) {
    if (num_elements < 0) {
        throw invalid_argument("Number of elements must be non-negative.");
    }
    Vector<T> vector(num_elements);
    for (int i = 0; i < num_elements; i++) {
        vector(i) = value;
    }
    return vector;
}

/**
 * @brief Create a vector of a given size with random values between [0, 1]
 * @param num_elements The number of elements in the vector.
 * @return The vector of random values.
 */
Vector<double> rand(int num_elements) {
    if (num_elements < 0) {
        throw invalid_argument("Number of elements must be non-negative.");
    }

    // Seed random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(0.0, 1.0);

    Vector<double> vector(num_elements);
    for (int i = 0; i < num_elements; i++) {
        vector(i) = dis(gen);
    }
    return vector;
}

/**
 * @brief Find the square root of each element in the vector.
 * @param vector The vector to take the square root of.
 * @return The vector with the square root of each element.
 */
template <typename T>
Vector<T> sqrt(const Vector<T>& vector) {
    Vector<T> result(vector.size());
    for (int i = 0; i < vector.size(); i++) {
        result(i) = sqrt(vector(i));
    }
    return result;
}

/**
 * @brief Accumulate the sum of all elements in the vector.
 * @return The sum of all elements in the vector.
 */
template <typename T>
T Vector<T>::accumulate() const {
    T sum = 0;
    for (int i = 0; i < n_elements; i++) {
        sum += data[i];
    }
    return sum;
}

/**
 * @brief Return the absolute value of each element in the vector.
 * @return The vector with the absolute value of each element.
*/
template <typename T>
Vector<T> abs(const Vector<T>& vector) {
    Vector<T> result(vector.size());
    for (int i = 0; i < vector.size(); i++) {
        result(i) = abs(vector(i));
    }
    return result;
}

/**
 * @brief Get the maximum element in the vector.
 * @return The maximum element in the vector.
*/
template <typename T>
T Vector<T>::max() const {
    T max = data[0];
    for (int i = 1; i < n_elements; i++) {
        if (data[i] > max) {
            max = data[i];
        }
    }
    return max;
}

// ------------------------ Template Instantiations ------------------------ //
template class Vector<int>;
template class Vector<float>;
template class Vector<double>;

template Vector<int> zeros(int num_elements);
template Vector<float> zeros(int num_elements);
template Vector<double> zeros(int num_elements);

template Vector<int> ones(int num_elements);
template Vector<float> ones(int num_elements);
template Vector<double> ones(int num_elements);

template Vector<int> fill(int num_elements, int value);
template Vector<float> fill(int num_elements, float value);
template Vector<double> fill(int num_elements, double value);

template Vector<int> sqrt(const Vector<int>& vector);
template Vector<float> sqrt(const Vector<float>& vector);
template Vector<double> sqrt(const Vector<double>& vector);

template Vector<int> abs(const Vector<int>& vector);
template Vector<float> abs(const Vector<float>& vector);
template Vector<double> abs(const Vector<double>& vector);