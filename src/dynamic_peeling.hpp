#ifndef FAST_MATRIX_MULTIPLICATION_DYNAMIC_PEELING_HPP
#define FAST_MATRIX_MULTIPLICATION_DYNAMIC_PEELING_HPP

#include <cassert>
#include "matrix.hpp"
#include "multiply_classic.hpp"

// performs dynamic peeling for matrix product <n, k, m> algorithm
// (dimensions n, k, and m mean in how many blocks we split matrices A and B)
// this function directly adds needed products to matrix C
template<typename Scalar>
void dynamic_peeling(const Matrix<Scalar> &A, const Matrix<Scalar> &B, Matrix<Scalar> &C,
                     unsigned int n, unsigned int k, unsigned int m) {
    // check if matrix dimensions are valid
    assert(A.rows == C.rows && A.cols == B.rows && B.cols == C.cols);

    // how many rows and cols were included in algorithm
    unsigned int included_rows_A = (A.rows / n) * n,
            included_cols_A = (A.cols / k) * k,
            included_rows_B = (B.rows / k) * k,
            included_cols_B = (B.cols / m) * m;

    // for how many rows do we need to do dynamic peeling?
    unsigned int need_peeling_rows_A = A.rows - included_rows_A,
            need_peeling_cols_A = A.cols - included_cols_A,
            need_peeling_rows_B = B.rows - included_rows_B,
            need_peeling_cols_B = B.cols - included_cols_B;

    // add product of not included columns from A * not included rows from B
    // we will split this in 3 blocks we can easily calculate

    // first block
    // A * B = C
    // | . . O |   | . . . |   | O O . |
    // | . . O | * | . . . | = | O O . |
    // | . . . |   | O O . |   | . . . |
    if (need_peeling_cols_A > 0) {
        Matrix<Scalar> A_extra = A.subblock({0, included_cols_A}, {included_rows_A, need_peeling_cols_A});
        Matrix<Scalar> B_extra = B.subblock({included_rows_B, 0}, {need_peeling_rows_B, included_cols_B});

        // calculate this product and add it to large matrix C in correct place
        Matrix<Scalar> product = multiply_classic(A_extra, B_extra);

        C.block_add({0, 0}, product);
    }

    // second block
    // A * B = C
    // | O O O |   | . . O |   | . . O |
    // | O O O | * | . . O | = | . . O |
    // | O O O |   | . . O |   | . . O |
    if (need_peeling_cols_B > 0) {
        Matrix<Scalar> B_extra = B.subblock({0, included_cols_B}, {B.rows, need_peeling_cols_B});

        // calculate product and add it to product C
        Matrix<Scalar> product = multiply_classic(A, B_extra);
        C.block_add({0, included_cols_B}, product);
    }

    // third block
    // A * B = C
    // | . . . |   | O O . |   | . . . |
    // | . . . | * | O O . | = | . . . |
    // | O O O |   | O O . |   | O O . |
    if (need_peeling_rows_A > 0) {
        Matrix<Scalar> A_extra = A.subblock({included_rows_A, 0}, {need_peeling_rows_A, A.cols});
        Matrix<Scalar> B_extra = B.subblock({0, 0}, {B.rows, included_cols_B});

        // calculate this product and add it to large matrix C in correct place
        Matrix<Scalar> product = multiply_classic(A_extra, B_extra);
        C.block_add({included_rows_A, 0}, product);
    }
}

#endif //FAST_MATRIX_MULTIPLICATION_DYNAMIC_PEELING_HPP
