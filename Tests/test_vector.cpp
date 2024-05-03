/**
 * @file test_vector.cpp
 * @brief Test initialization of the Vector Class.
 * 
 * This file contains the test cases for creating the Vector class where
 * normal and edge cases are tested to ensure proper function.
 * 
 * @author Joshua Blomgren
 * @date December 8, 2023
 */

#include "vector.h"
#include <cassert>

/**
 * @brief Test the default constructor.
 */
void test_default_constructor() {
    cout << "Testing the default constructor..." << endl;
    Vector<int> vector;
    assert(vector.size() == 0);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the constructor with given size.
 */
void test_constructor_with_size() {
    cout << "Testing the constructor with given size..." << endl;
    Vector<int> vector(10);
    assert(vector.size() == 10);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the copy constructor.
 */
void test_copy_constructor() {
    cout << "Testing the copy constructor..." << endl;
    Vector<int> vector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Vector<int> vector_copy(vector);
    assert(vector_copy == vector);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the () operator.
 */
void test_assignment_operator() {
    cout << "Testing the () operator..." << endl;
    Vector<int> vector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    assert(vector(0) == 1);
    assert(vector(1) == 2);
    assert(vector(2) == 3);
    assert(vector(3) == 4);
    assert(vector(4) == 5);
    assert(vector(5) == 6);
    assert(vector(6) == 7);
    assert(vector(7) == 8);
    assert(vector(8) == 9);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the size() function.
 */
void test_size() {
    cout << "Testing the size() function..." << endl;
    Vector<int> vector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    assert(vector.size() == 9);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test printing the vector.
 */
void test_print() {
    cout << "Testing printing the vector..." << endl;
    Vector<int> vector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    vector.print();
    cout << "Test passed!" << endl;
}

/**
 * @brief Test resizing the vector.
 */
void test_resize() {
    cout << "Testing resizing the vector..." << endl;
    Vector<int> vector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    vector.resize(5);
    assert(vector.size() == 5);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test resetting the size of the vector.
 */
void test_reset_size() {
    cout << "Testing resetting the size of the vector..." << endl;
    Vector<int> vector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    vector.reset_size(5);
    assert(vector.size() == 5);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test overloading the << operator to print the vector.
 */
void test_overload_print() {
    cout << "Testing overloading the << operator to print the vector..." << endl;
    Vector<int> vector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    cout << vector;
    cout << "Test passed!" << endl;
}

/**
 * @brief Test overloading the = operator to assign one vector to another.
 */
void test_overload_assignment() {
    cout << "Testing overloading the = operator to assign one vector to another..." << endl;
    Vector<int> vector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Vector<int> vector2 = vector;
    assert(vector2 == vector);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test overloading the == operator to check if two vectors are equal.
 */
void test_overload_equality() {
    cout << "Testing overloading the == operator to check if two vectors are equal..." << endl;
    Vector<int> vector = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    Vector<int> vector2 = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    assert(vector == vector2);
    cout << "Test passed!" << endl;
}

/**
 * @brief Test the fill functions (zeros, ones, fill).
 */
void test_fill() {
    cout << "Testing the zeros function..." << endl;
    Vector<int> vector = zeros<int>(10);
    assert(vector == Vector<int>({0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));
    
    cout << "Testing the ones function..." << endl;
    vector = ones<int>(10);
    assert(vector == Vector<int>({1, 1, 1, 1, 1, 1, 1, 1, 1, 1}));

    cout << "Testing the fill function..." << endl;
    Vector<int> five_vector = fill(10, 5);
    assert(five_vector == Vector<int>({5, 5, 5, 5, 5, 5, 5, 5, 5, 5}));
    cout << "Test passed!" << endl;
}

/**
 * @brief Test random vector generation.
 */
void test_rand() {
    cout << "Testing random vector generation..." << endl;
    Vector<double> vector = rand(5);
    cout << vector;
    cout << "Test passed!" << endl;
}

int main() {
    test_default_constructor();
    test_constructor_with_size();
    test_copy_constructor();
    test_assignment_operator();
    test_size();
    test_resize();
    test_reset_size();
    test_print();
    test_overload_print();
    test_overload_assignment();
    test_overload_equality();
    test_fill();
    test_rand();

    return 0;
}