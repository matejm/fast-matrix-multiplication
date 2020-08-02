#ifndef FAST_MATRIX_MULTIPLICATION_MULTIPLY_CLASSIC_HPP
#define FAST_MATRIX_MULTIPLICATION_MULTIPLY_CLASSIC_HPP

#include <cassert>
#include <vector>
#include "matrix.hpp"

template<class Scalar>
Matrix<Scalar> multiply_classic(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    // check dimensions
    assert(A.cols == B.rows);

    // create new matrix
    Matrix<Scalar> C = Matrix<Scalar>::zeros(A.rows, B.cols);

    for (unsigned int i = 0; i < A.rows; ++i) {
        for (unsigned int k = 0; k < A.cols; ++k) {
            // calculate c_ij
            auto Aik = A[{i, k}];
            // this will be Bkj, but using iterators improves multiplication speed by a factor of 2
            auto iter_B = B.data.begin() + k * B.cols;
            // this will be Cij
            auto iter_C = C.data.begin() + i * C.cols;
            for (unsigned int j = 0; j < B.cols; ++j) {
                *(iter_C + j) += Aik * (*(iter_B + j));
            }
        }
    }
    return C;
}

#endif //FAST_MATRIX_MULTIPLICATION_MULTIPLY_CLASSIC_HPP
