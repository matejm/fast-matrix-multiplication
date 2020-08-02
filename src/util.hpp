#ifndef FAST_MATRIX_MULTIPLICATION_UTIL_HPP
#define FAST_MATRIX_MULTIPLICATION_UTIL_HPP

#include <iostream>
#include <random>
#include <string>
#include "matrix.hpp"

std::default_random_engine generator;

// Simple function that prints a matrix.
template<class Scalar>
void print_matrix(const Matrix<Scalar> &A, const std::string &matrix_name = "") {
    if (matrix_name != "") {
        std::cout << matrix_name << " ";
    } else {
        std::cout << "matrix ";
    }
    std::cout << A;
}

template<typename Scalar>
void print_differences(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    assert(A.rows == B.rows && A.cols == B.cols);
    std::vector<bool> data(A.data.size());

    for (unsigned int i = 0; i < A.rows; i++) {
        for (unsigned int j = 0; j < A.cols; j++) {
            data[i * A.cols + j] = A[{i, j}] == B[{i, j}];
        }
    }

    print_matrix(Matrix<bool>(data, A.rows, A.cols));
}

// constructs matrix of size rows x cols, filled with random int values,
// values are from 0 to max_size (inclusive)
Matrix<int> random_int_matrix(unsigned int rows, unsigned int cols, unsigned int max_size = 10) {
    // create new random generator
    std::uniform_int_distribution<int> distribution(0, max_size);

    // create new matrix data and fill it
    std::vector<int> data(rows * cols);

    for (unsigned int i = 0; i < rows * cols; ++i) {
        data[i] = distribution(generator);
    }

    return Matrix<int>(data, rows, cols);
}

// constructs matrix of size rows x cols, filled with random int values,
// values are from 0 to max_size (inclusive)
Matrix<double> random_float_matrix(unsigned int rows, unsigned int cols, unsigned int max_size = 10) {
    // create new random generator
    std::uniform_int_distribution<int> distribution(0, max_size);

    // create new matrix data and fill it
    std::vector<double> data(rows * cols);

    for (unsigned int i = 0; i < rows * cols; ++i) {
        data[i] = distribution(generator);
    }

    return Matrix<double>(data, rows, cols);
}


#endif //FAST_MATRIX_MULTIPLICATION_UTIL_HPP
