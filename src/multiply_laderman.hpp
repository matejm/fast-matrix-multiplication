#ifndef FAST_MATRIX_MULTIPLICATION_MULTIPLY_LADERMAN_HPP
#define FAST_MATRIX_MULTIPLICATION_MULTIPLY_LADERMAN_HPP

#include "matrix.hpp"
#include "multiply_classic.hpp"
#include "dynamic_peeling.hpp"

// when to switch to classic matrix multiplication
// for Laderman's algorithm, this has to be at least 2
const unsigned int laderman_threshold = 200;

// Julian B. Larderman: A Noncomutative Algorithm for Multiplying 3x3 Matrices Using 23 Multiplications
// similar to Strassen, except it works on 3x3 matrices.
// static and dynamic versions could be implemented, but adding zeros to next power of 3
// is not very practical, this is why only dynamic peeling version was implemented
template<typename Scalar>
Matrix<Scalar> multiply_laderman(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    // dimension check
    assert(A.cols == B.rows);

    // if any of the dimensions is too small, laderman's algorithm wont help
    if (std::min(A.rows, std::min(A.cols, B.cols)) <= laderman_threshold) {
        return multiply_classic(A, B);
    }

    // subblock sizes
    std::pair<unsigned int, unsigned int>
            block_A = {A.rows / 3, A.cols / 3},
            block_B = {B.rows / 3, B.cols / 3},
            product_block = {A.rows / 3, B.cols / 3};
    // block dimensions
    unsigned int block_rows_A = block_A.first,
            block_cols_A = block_A.second,
            block_rows_B = block_B.first,
            block_cols_B = block_B.second,
            product_rows = product_block.first,
            product_cols = product_block.second;

    // split matrices in subblocks
    Matrix<Scalar> A11 = A.subblock({0, 0}, block_A),
            A12 = A.subblock({0, block_cols_A}, block_A),
            A13 = A.subblock({0, 2 * block_cols_A}, block_A),
            A21 = A.subblock({block_rows_A, 0}, block_A),
            A22 = A.subblock({block_rows_A, block_cols_A}, block_A),
            A23 = A.subblock({block_rows_A, 2 * block_cols_A}, block_A),
            A31 = A.subblock({2 * block_rows_A, 0}, block_A),
            A32 = A.subblock({2 * block_rows_A, block_cols_A}, block_A),
            A33 = A.subblock({2 * block_rows_A, 2 * block_cols_A}, block_A);

    Matrix<Scalar> B11 = B.subblock({0, 0}, block_B),
            B12 = B.subblock({0, block_cols_B}, block_B),
            B13 = B.subblock({0, 2 * block_cols_B}, block_B),
            B21 = B.subblock({block_rows_B, 0}, block_B),
            B22 = B.subblock({block_rows_B, block_cols_B}, block_B),
            B23 = B.subblock({block_rows_B, 2 * block_cols_B}, block_B),
            B31 = B.subblock({2 * block_rows_B, 0}, block_B),
            B32 = B.subblock({2 * block_rows_B, block_cols_B}, block_B),
            B33 = B.subblock({2 * block_rows_B, 2 * block_cols_B}, block_B);

    // create larger matrix for result
    Matrix<Scalar> C = Matrix<Scalar>::zeros(A.rows, B.cols);

    // temporary matrix, here we will store products
    Matrix<Scalar> P;

    // 23 products

    // P1 = (A11 + A12 + A13 - A21 - A22 - A32) B22
    P = multiply_laderman(A11 + A12 + A13 - A21 - A22 - A32 - A33, B22);
    // add P1 to C12
    C.block_add({0, product_cols}, P);

    // P2 = (A11 - A21) (B22 - B12)
    P = multiply_laderman(A11 - A21, B22 - B12);
    // add P2 to C21, C22
    C.block_add({product_rows, 0}, P);
    C.block_add({product_rows, product_cols}, P);

    // P3 = A22 (B12 - B11 + B21 - B22 - B23 - B31 + B33)
    P = multiply_laderman(A22, B12 - B11 + B21 - B22 - B23 - B31 + B33);
    // add P3 to C21
    C.block_add({product_rows, 0}, P);

    // P4 = (A21 - A11 + A22) (B11 - B12 + B22)
    P = multiply_laderman(A21 - A11 + A22, B11 - B12 + B22);
    // add P4 to C12, C21, C22
    C.block_add({0, product_cols}, P);
    C.block_add({product_rows, 0}, P);
    C.block_add({product_rows, product_cols}, P);

    // P5 = (A21 + A22) (B12 - B11)
    P = multiply_laderman(A21 + A22, B12 - B11);
    // add P5 to C12, C22
    C.block_add({0, product_cols}, P);
    C.block_add({product_rows, product_cols}, P);

    // P6 = A11 B11
    P = multiply_laderman(A11, B11);
    // add P6 to C11, C12, C13, C21, C22, C31, C33
    C.block_add({0, 0}, P);
    C.block_add({0, product_cols}, P);
    C.block_add({0, 2 * product_cols}, P);
    C.block_add({product_rows, 0}, P);
    C.block_add({product_rows, product_cols}, P);
    C.block_add({2 * product_rows, 0}, P);
    C.block_add({2 * product_rows, 2 * product_cols}, P);

    // P7 = (A31 - A11 + A32) (B11 - B13 + B23)
    P = multiply_laderman(A31 - A11 + A32, B11 - B13 + B23);
    // add P7 to C13, C31, C33
    C.block_add({0, 2 * product_cols}, P);
    C.block_add({2 * product_rows, 0}, P);
    C.block_add({2 * product_rows, 2 * product_cols}, P);

    // P8 = (A31 - A11) (B13 - B23)
    P = multiply_laderman(A31 - A11, B13 - B23);
    // add P8 to C31, C33
    C.block_add({2 * product_rows, 0}, P);
    C.block_add({2 * product_rows, 2 * product_cols}, P);

    // P9 = (A31 + A32) (B13 - B23)
    P = multiply_laderman(A31 + A32, B13 - B11);
    // add P9 to C13, C33
    C.block_add({0, 2 * product_cols}, P);
    C.block_add({2 * product_rows, 2 * product_cols}, P);

    // P10 = (A11 + A12 + A13 - A22 - A23 - A31 - A32) B23
    P = multiply_laderman(A11 + A12 + A13 - A22 - A23 - A31 - A32, B23);
    // add P10 to C13
    C.block_add({0, 2 * product_cols}, P);

    // P11 = A32 (B13 - B11 + B21 - B22 - B23 - B31 + B32)
    P = multiply_laderman(A32, B13 - B11 + B21 - B22 - B23 - B31 + B32);
    // add P11 to C31
    C.block_add({2 * product_rows, 0}, P);

    // P12 = (A32 - A13 + A33) (B22 + B31 - B32)
    P = multiply_laderman(A32 - A13 + A33, B22 + B31 - B32);
    // add P12 to C12, C31, C32
    C.block_add({0, product_cols}, P);
    C.block_add({2 * product_rows, 0}, P);
    C.block_add({2 * product_rows, product_cols}, P);

    // P13 = (A13 - A33) (B22 - B23)
    P = multiply_laderman(A13 - A33, B22 - B32);
    // add P13 to C31, C32
    C.block_add({2 * product_rows, 0}, P);
    C.block_add({2 * product_rows, product_cols}, P);

    // P14 = A13 B31
    P = multiply_laderman(A13, B31);
    // add P14 to C11, C12, C13, C21, C23, C31, C32
    C.block_add({0, 0}, P);
    C.block_add({0, product_cols}, P);
    C.block_add({0, 2 * product_cols}, P);
    C.block_add({product_rows, 0}, P);
    C.block_add({product_rows, 2 * product_cols}, P);
    C.block_add({2 * product_rows, 0}, P);
    C.block_add({2 * product_rows, product_cols}, P);

    // P15 = (A32 + A33) (B32 - B31)
    P = multiply_laderman(A32 + A33, B32 - B31);
    // add P15 to C12, C32
    C.block_add({0, product_cols}, P);
    C.block_add({2 * product_rows, product_cols}, P);

    // P16 = (A22 - A13 + A23) (B23 + B31 - B33)
    P = multiply_laderman(A22 - A13 + A23, B23 + B31 - B33);
    // add P16 to C13, C21, C23
    C.block_add({0, 2 * product_cols}, P);
    C.block_add({product_rows, 0}, P);
    C.block_add({product_rows, 2 * product_cols}, P);

    // P17 = (A13 - A23) (B23 - B33)
    P = multiply_laderman(A13 - A23, B23 - B33);
    // add P17 to C21, C23
    C.block_add({product_rows, 0}, P);
    C.block_add({product_rows, 2 * product_cols}, P);

    // P18 = (A22 + A23) (B33 - B31)
    P = multiply_laderman(A22 + A23, B33 - B31);
    // add P18 to C13, C23
    C.block_add({0, 2 * product_cols}, P);
    C.block_add({product_rows, 2 * product_cols}, P);

    // P19 = A12 B21
    P = multiply_laderman(A12, B21);
    // add P19 to C11
    C.block_add({0, 0}, P);

    // P20 = A23 B32
    P = multiply_laderman(A23, B32);
    // add P20 to C22
    C.block_add({product_rows, product_cols}, P);

    // P21 = A21 B13
    P = multiply_laderman(A21, B13);
    // add P21 to C23
    C.block_add({product_rows, 2 * product_cols}, P);

    // P22 = A31 B12
    P = multiply_laderman(A31, B12);
    // add P22 to C32
    C.block_add({2 * product_rows, product_cols}, P);

    // P23 = A33 B33
    P = multiply_laderman(A33, B33);
    // add P23 to C33
    C.block_add({2 * product_rows, 2 * product_cols}, P);

    // fix remaining row and column if dimensions are odd
    dynamic_peeling(A, B, C, 3, 3, 3);

    return C;

}

#endif //FAST_MATRIX_MULTIPLICATION_MULTIPLY_LADERMAN_HPP
