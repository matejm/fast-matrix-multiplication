#include "helpers.hpp"

// constructs matrix of size rows x cols, filled with random int values,
// values are from 0 to max_size (inclusive)
Matrix<int> random_int_matrix(unsigned int rows, unsigned int cols, unsigned int max_size) {
    // create new random generator
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(0, max_size);

    // create new matrix data and fill it
    std::vector<int> data(rows * cols);

    for (unsigned int i = 0; i < rows * cols; ++i) {
        data[i] = distribution(generator);
    }

    return Matrix<int>(data, rows, cols);
}

// constructs matrix of size rows x cols, filled with random int values,
// values are from 1 to max_size (inclusive)
Matrix<double> random_float_matrix(unsigned int rows, unsigned int cols, unsigned int max_size) {
    // create new random generator
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1, max_size);

    // create new matrix data and fill it
    std::vector<double> data(rows * cols);

    for (unsigned int i = 0; i < rows * cols; ++i) {
        data[i] = distribution(generator);
    }

    return Matrix<double>(data, rows, cols);
}

// one of possible measures how far approximate matrix is from correct matrix
double maximum_relative_difference(const Matrix<double> &correct, const Matrix<double> &approx) {
    // check dimensions
    assert(correct.rows == approx.rows && correct.cols == approx.cols);

    // vector of differences
    std::vector<double> differences = (correct - approx).data;

    double max_value = 0;

    // calculate relative differences, divide by correct value
    for (unsigned int i = 0; i < differences.size(); ++i) {
        max_value = std::max(
                max_value,
                std::abs(differences[i] / correct.data[i])
        );
    }

    return max_value;
}
