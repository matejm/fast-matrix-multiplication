#ifndef FAST_MATRIX_MULTIPLICATION_MULTIPLY_BINI_HPP
#define FAST_MATRIX_MULTIPLICATION_MULTIPLY_BINI_HPP

#include <cassert>
#include "matrix.hpp"
#include "dynamic_peeling.hpp"
#include "polynomial.hpp"

const unsigned int bini_threshold = 200;

// Exact multiplication, problem is that multiplication of polynomials is not O(1) anymore.
template<typename Scalar>
Matrix<Scalar> multiply_bini_exact(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    // convert to polynomials
    Matrix<Polynomial<Scalar>> poly_A(A);
    Matrix<Polynomial<Scalar>> poly_B(B);

    Matrix<Polynomial<Scalar>> poly_C = multiply_bini(poly_A, poly_B, Polynomial<Scalar>::epsilon());

    // convert back from polynomials
    return polynomial_to_scalar(poly_C);
}

// actual Bini's algorithm
template<typename Poly>
Matrix<Poly> multiply_bini(const Matrix<Poly> &A, const Matrix<Poly> &B, const Poly &epsilon) {
    // check dimensions
    assert(A.cols == B.rows);

    // if matrices are too small for Bini's algorithm we have nothing to do but multiply it classicaly
    if (A.rows < 2 || A.cols < 2 || B.cols < 3) {
        return multiply_classic(A, B);
    } else if (A.rows <= bini_threshold || A.cols <= bini_threshold || B.cols <= bini_threshold) {
        return multiply_classic(A, B);
    }

    // create subblocks
    // last line and up to 2 last columns may not be included, this is handled by dynamic peeling
    unsigned int block_rows_A = A.rows / 2,
            block_cols_A = A.cols / 2,
            block_rows_B = B.rows / 2,
            block_cols_B = B.cols / 3;

    std::pair<unsigned int, unsigned int>
            block_A = {block_rows_A, block_cols_A},
            block_B = {block_rows_B, block_cols_B};

    // | A11 A12 | | B11 B12 B13 |
    // | A21 A22 | | B21 B22 B23 |
    Matrix<Poly> A11 = A.subblock({0, 0}, block_A),
            A12 = A.subblock({0, block_cols_A}, block_A),
            A21 = A.subblock({block_rows_A, 0}, block_A),
            A22 = A.subblock({block_rows_A, block_cols_A}, block_A);

    Matrix<Poly> B11 = B.subblock({0, 0}, block_B),
            B12 = B.subblock({0, block_cols_B}, block_B),
            B13 = B.subblock({0, 2 * block_cols_B}, block_B),
            B21 = B.subblock({block_rows_B, 0}, block_B),
            B22 = B.subblock({block_rows_B, block_cols_B}, block_B),
            B23 = B.subblock({block_rows_B, 2 * block_cols_B}, block_B);

    // create new empty matrix for product
    // | C11 C12 C13 |
    // | C21 C22 C23 |
    Matrix<Poly> C = Matrix<Poly>::zeros(A.rows, B.cols);

    // dimensions of a block in a product
    unsigned int block_rows_C = block_rows_A, block_cols_C = block_cols_B;

    // matrix for storing products
    Matrix<Poly> P;

    // first half
    // | C11 C12 |
    // | C21 ... |

    // P1 = (A12 + e A22) B21
    P = multiply_bini(A12 + epsilon * A22, B21, epsilon);
    // Add e P1 to C11 and P1 to C21
    C.block_add({0, 0}, epsilon * P);
    C.block_add({block_rows_C, 0}, P);

    // P2 = A11 (B11 + e B12)
    P = multiply_bini(A11, B11 + epsilon * B12, epsilon);
    // Add e P2 to C11 and P2 to C12
    C.block_add({0, 0}, epsilon * P);
    C.block_add({0, block_cols_C}, P);

    // P3 = A12 (B11 + B21 + e B22)
    P = multiply_bini(A12, B11 + B21 + epsilon * B22, epsilon);
    // Subtract P3 from C21
    C.block_subtract({block_rows_C, 0}, P);

    // P4 = (A11 + A12 + e A21) B11
    P = multiply_bini(A11 + A12 + epsilon * A21, B11, epsilon);
    // Subtract P4 from C12
    C.block_subtract({0, block_cols_C}, P);

    // P5 = (A12 + e A21) (B11 + e B22)
    P = multiply_bini(A12 + epsilon * A21, B11 + epsilon * B22, epsilon);
    // Add P5 to C12 and C21
    C.block_add({0, block_cols_C}, P);
    C.block_add({block_rows_C, 0}, P);

    // second half
    // | ... C13 |
    // | C22 C23 |
    // how to get to formulas for P1, ... P5: transpose upper matrix and rename indices.

    // transpose all blocks (we won't need B11 and B21, we do not need to transpose them).
    A11 = A11.transposed();
    A12 = A12.transposed();
    A21 = A21.transposed();
    A22 = A22.transposed();
    B12 = B12.transposed();
    B22 = B22.transposed();
    B13 = B13.transposed();
    B23 = B23.transposed();

    // P1 = (B13 + e B12) A21
    // After multiplying, transposed back
    P = multiply_bini(B13 + epsilon * B12, A21, epsilon).transposed();
    // Add e P1 to C23 and P1 to C22
    C.block_add({block_rows_C, 2 * block_cols_C}, epsilon * P);
    C.block_add({block_rows_C, block_cols_C}, P);

    // P2 = B23 (A22 + e A12)
    P = multiply_bini(B23, A22 + epsilon * A12, epsilon).transposed();
    // Add e P2 to C23 and P2 to C13
    C.block_add({block_rows_C, 2 * block_cols_C}, epsilon * P);
    C.block_add({0, 2 * block_cols_C}, P);

    // P3 = B13 (A22 + A21 + e A11)
    P = multiply_bini(B13, A22 + A21 + epsilon * A11, epsilon).transposed();
    // Subtract P3 from C22
    C.block_subtract({block_rows_C, block_cols_C}, P);

    // P4 = (B23 + B13 + e B22) A22
    P = multiply_bini(B23 + B13 + epsilon * B22, A22, epsilon).transposed();
    // Subtract P4 from C13
    C.block_subtract({0, 2 * block_cols_C}, P);

    // P5 = (B13 + e B22) (A22 + e A11)
    P = multiply_bini(B13 + epsilon * B22, A22 + epsilon * A11, epsilon).transposed();
    // Add P5 to C13 and C22
    C.block_add({0, 2 * block_cols_C}, P);
    C.block_add({block_rows_C, block_cols_C}, P);

    // using bini's algorithm now we got epsilon * C, now we have to divide by epsilon.
    C /= epsilon;

    // dynamic peeling for not included rows and cols
    dynamic_peeling(A, B, C, 2, 2, 3);

    return C;
}


#endif //FAST_MATRIX_MULTIPLICATION_MULTIPLY_BINI_HPP
