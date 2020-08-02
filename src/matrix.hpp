#ifndef FAST_MATRIX_MULTIPLICATION_MATRIX_HPP
#define FAST_MATRIX_MULTIPLICATION_MATRIX_HPP

#include<iostream>
#include<array>
#include<vector>
#include<cassert>
#include<algorithm>
#include<iterator>
#include <tuple>

template<class Scalar>
class Matrix {
public:
    // number of rows in this matrix
    unsigned int rows;
    // number of columns in this matrix
    unsigned int cols;

    // where actual matrix data is stored
    // (maybe this should be changed to C's array which is faster)
    std::vector<Scalar> data;

    // default constructor is empty matrix
    Matrix() : rows(0), cols(0), data({}) {}

    // constructs matrix of size rows x cols, filled with initial_value
    Matrix(const unsigned int _rows, const unsigned int _cols, const Scalar &initial_value) : rows(_rows), cols(_cols) {
        // allocate necessary space and fill
        data = std::vector<Scalar>(_rows * _cols, initial_value);
    }

    // constructs matrix of size rows x cols, filled with zeros
    static const Matrix<Scalar> zeros(const unsigned int rows, const unsigned int cols) {
        return Matrix(rows, cols, Scalar(0));
    }

    // copy matrix from 1D array
    // no dimension checks are performed
    Matrix(const std::vector<Scalar> &_data, const unsigned int _rows, const unsigned int _cols) : data(_data),
                                                                                                   rows(_rows),
                                                                                                   cols(_cols) {}

    // Cast matrix of other type to this type.
    template<typename OtherScalar>
    Matrix(const Matrix<OtherScalar> A) {
        rows = A.rows;
        cols = A.cols;
        data = std::vector<Scalar>(A.data.begin(), A.data.end());
    }

    // matrix is indexed from (0, 0) to (rows-1, cols-1)
    // because reference is returned, values can also be set using this
    Scalar &operator[](const std::pair<unsigned int, unsigned int> location) {
        const int i = location.first;
        const int j = location.second;

        return data[i * cols + j];
    }

    // same as above, except this is const version
    const Scalar &operator[](const std::pair<unsigned int, unsigned int> location) const {
        const int i = location.first;
        const int j = location.second;

        return data[i * cols + j];
    }

    // return transposed  matrix
    Matrix<Scalar> transposed() {
        std::vector<Scalar> new_data(rows * cols);
        for (unsigned int i = 0; i < rows; ++i) {
            for (unsigned int j = 0; j < cols; ++j) {
                new_data[j * rows + i] = data[i * cols + j];
            }
        }
        return Matrix<Scalar>(new_data, cols, rows);
    }

    // block sub-matrix of current matrix
    // create new smaller matrix and copy data to it
    Matrix<Scalar> subblock(std::pair<unsigned int, unsigned int> top_left,
                            std::pair<unsigned int, unsigned int> block_size) const {
        // unpack
        unsigned int start_row, start_col;
        unsigned int block_rows, block_cols;

        std::tie(start_row, start_col) = top_left;
        std::tie(block_rows, block_cols) = block_size;

        // check dimensions
        assert(start_row + block_rows <= rows);
        assert(start_col + block_cols <= cols);

        // allocate new vector and create new matrix
        std::vector<Scalar> block_data(block_rows * block_cols);
        Matrix<Scalar> block(block_data, block_rows, block_cols);

        // get iterator of newly created matrix
        auto output_iter = block.data.begin();

        for (unsigned int i = 0; i < block_rows; ++i) {
            // copy i-th row to newly created
            auto input_iter = data.begin() + ((start_row + i) * cols + start_col);
            for (unsigned int j = 0; j < block_cols; ++j) {
                // move value
                *output_iter++ += *input_iter++;
            }
        }
        return block;
    }

    // add block starting from top_left to this matrix
    Matrix<Scalar> &block_add(std::pair<unsigned int, unsigned int> top_left, const Matrix<Scalar> &block) {
        unsigned int start_row = top_left.first;
        unsigned int start_col = top_left.second;

        // check if dimensions are correct
        assert(start_row + block.rows <= rows && start_col + block.cols <= cols);

        // get iterator of block matrix
        auto block_iter = block.data.begin();

        for (unsigned int i = 0; i < block.rows; ++i) {
            // add i-th row to our data
            auto data_iter = data.begin() + ((start_row + i) * cols + start_col);
            for (unsigned int j = 0; j < block.cols; ++j) {
                // add value
                *data_iter++ += *block_iter++;
            }
        }
        return *this;
    }

    // subtract block starting from top_left to this matrix
    // exactly same as block_add except here we are subtracting
    Matrix<Scalar> &block_subtract(std::pair<unsigned int, unsigned int> top_left, const Matrix<Scalar> &block) {
        unsigned int start_row = top_left.first;
        unsigned int start_col = top_left.second;

        // check if dimensions are correct
        assert(start_row + block.rows <= rows && start_col + block.cols <= cols);

        // get iterator of block matrix
        auto block_iter = block.data.begin();

        for (unsigned int i = 0; i < block.rows; ++i) {
            // add i-th row to our data
            auto data_iter = data.begin() + ((start_row + i) * cols + start_col);
            for (unsigned int j = 0; j < block.cols; ++j) {
                // subtract value
                *data_iter++ -= *block_iter++;
            }
        }
        return *this;
    }

    // add other matrix to this matrix
    // does not create new matrix, changes current one
    Matrix<Scalar> &operator+=(const Matrix<Scalar> &other) {
        // check if dimensions are correct
        assert(rows == other.rows && cols == other.cols);

        // add other's data to our data
        auto iter = other.data.begin();
        for (auto &element : data) {
            element += *iter;
            iter++;
        }
        return *this;
    }

    // subtract other matrix from this matrix
    // does not create new matrix, changes current one
    Matrix<Scalar> &operator-=(const Matrix<Scalar> &other) {
        // ensure that both matrices have same shape
        assert(rows == other.rows && cols == other.cols);

        // subtract other's data from our data
        auto iter = other.data.begin();
        for (auto &element : data) {
            element -= *iter;
            iter++;
        }
        return *this;
    }

    // multiply matrix with a scalar
    Matrix<Scalar> &operator*=(const Scalar &s) {
        // multiply all elements with a scalar value
        for (auto &element : data) {
            element *= s;
        }
        return *this;
    }

    // divide matrix by a scalar
    Matrix<Scalar> &operator/=(const Scalar &s) {
        // multiply all elements with a scalar value
        for (auto &element : data) {
            element /= s;
        }
        return *this;
    }

    // check equality by elements
    bool operator==(const Matrix<Scalar> &other) const {
        if (rows != other.rows || cols != other.cols) return false;

        // if dimensions are correct, we need to check all items
        for (unsigned int i = 0; i < rows * cols; ++i) {
            if (data[i] != other.data[i]) return false;
        }
        return true;
    }

    bool operator!=(const Matrix<Scalar> &other) const {
        return !(*this == other);
    }
};

// Adds two matrices, creates new matrix object, leaves original matrices unchanged
template<class Scalar>
Matrix<Scalar> operator+(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    Matrix<Scalar> result = A;
    result += B;
    return result;
}

// Subtracts one matrix from another, creates new matrix object, leaves original matrices unchanged
template<class Scalar>
Matrix<Scalar> operator-(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    Matrix<Scalar> result = A;
    result -= B;
    return result;
}

// Multiply matrix with a scalar
template<class Scalar>
Matrix<Scalar> operator*(const Scalar &s, const Matrix<Scalar> &A) {
    Matrix<Scalar> result = A;
    result *= s;
    return result;
}

template<class Scalar>
Matrix<Scalar> operator==(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    return A == B;
}

// make matrix printable
template<typename Scalar>
std::ostream &operator<<(std::ostream &os, const Matrix<Scalar> &A) {
    auto iter = A.data.begin();

    os << "(" << A.rows << " x " << A.cols << ")" << std::endl;
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            os << *iter << '\t';
            iter++;
        }
        os << std::endl;
    }
    return os;
}

#endif //FAST_MATRIX_MULTIPLICATION_MATRIX_HPP
