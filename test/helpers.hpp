#ifndef FAST_MATRIX_MULTIPLICATION_HELPERS_HPP
#define FAST_MATRIX_MULTIPLICATION_HELPERS_HPP

#include <iostream>
#include <random>
#include <algorithm>
#include "matrix.hpp"

// constructs matrix of size rows x cols, filled with random int values,
// values are from 0 to max_size (inclusive)
Matrix<int> random_int_matrix(unsigned int rows, unsigned int cols, unsigned int max_size = 10);

// constructs matrix of size rows x cols, filled with random int values,
// values are from 0 to max_size (inclusive)
Matrix<double> random_float_matrix(unsigned int rows, unsigned int cols, unsigned int max_size = 10);

// one of possible measures how far approximate matrix is from correct matrix
double maximum_relative_difference(const Matrix<double> &correct, const Matrix<double> &approx);


#endif //FAST_MATRIX_MULTIPLICATION_HELPERS_HPP
