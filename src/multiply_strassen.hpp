#ifndef FAST_MATRIX_MULTIPLICATION_MULTIPLY_STRASSEN_HPP
#define FAST_MATRIX_MULTIPLICATION_MULTIPLY_STRASSEN_HPP

#include <algorithm>
#include "matrix.hpp"
#include "multiply_classic.hpp"
#include "dynamic_peeling.hpp"

const unsigned int strassen_threshold = 200;

// finds power of 2 larger (or same as) given value
unsigned int next_power_of_2(unsigned int value) {
    unsigned int power = 1;
    unsigned int n = 2;

    while (n < value) {
        power++;
        n *= 2;

        // otherwise we have too large matrices,
        // using integers for sizes is not sufficient
        assert(power < 32);
    }

    return n;
}

// strassen multiply using static padding:
// resize matrices A and B to the next power of 2
template<class Scalar>
Matrix<Scalar> multiply_strassen_static(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    // dimension check
    assert(A.cols == B.rows);

    // fill matrices A and B to both be squares with sizes powers of 2.
    unsigned int n = next_power_of_2(std::max(A.rows, std::max(B.cols, B.rows)));

    // create new large enough matrices
    Matrix<Scalar> new_A = Matrix<Scalar>::zeros(n, n);
    Matrix<Scalar> new_B = Matrix<Scalar>::zeros(n, n);

    // fill new matrix with values from A
    new_A.block_add({0, 0}, A);
    // fill new matrix with values from B
    new_B.block_add({0, 0}, B);

    // actual multiplication
    Matrix<Scalar> product = strassen(new_A, new_B);

    // return correct subblock (crop zeros)
    return product.subblock({0, 0}, {A.rows, B.cols});
}

// actual strassen multiplication
// here we can assume, that all matrices are square and all blocks are powers of 2,
// both matrices should have same size
template<class Scalar>
Matrix<Scalar> strassen(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    // (all dimensions should be the same)
    assert(A.rows == A.cols && A.cols == B.rows && B.rows == B.cols);

    const unsigned int size = A.rows;

    if (size <= strassen_threshold) {
        return multiply_classic(A, B);
    }

    // new block size, half of the original size
    unsigned int new_size = size / 2;
    std::pair<unsigned int, unsigned int> block = {new_size, new_size};

    // generate subblocks
    auto A11 = A.subblock({0, 0}, block);
    auto A12 = A.subblock({0, new_size}, block);
    auto A21 = A.subblock({new_size, 0}, block);
    auto A22 = A.subblock({new_size, new_size}, block);

    auto B11 = B.subblock({0, 0}, block);
    auto B12 = B.subblock({0, new_size}, block);
    auto B21 = B.subblock({new_size, 0}, block);
    auto B22 = B.subblock({new_size, new_size}, block);

    // create larger matrix for result
    Matrix<Scalar> C = Matrix<Scalar>::zeros(size, size);

    // temporary matrix, here we will store products
    Matrix<Scalar> P = Matrix<Scalar>::zeros(new_size, new_size);

    // P1 = (A11 + A22) (B11 + B22)
    P = strassen(A11 + A22, B11 + B22);
    // add P1 to C11 and C22
    C.block_add({0, 0}, P);
    C.block_add({new_size, new_size}, P);

    // P2 = (A21+ A22) B11
    P = strassen(A21 + A22, B11);
    // add P2 to C21 and subtract from C22
    C.block_add({new_size, 0}, P);
    C.block_subtract({new_size, new_size}, P);

    // P3 = A11 (B12 - B22)
    P = strassen(A11, B12 - B22);
    // add P3 to C12 and C22
    C.block_add({0, new_size}, P);
    C.block_add({new_size, new_size}, P);

    // P4 = A22 (-B11 + B21)
    P = strassen(A22, B21 - B11);
    // add P4 to C11 and C21
    C.block_add({0, 0}, P);
    C.block_add({new_size, 0}, P);

    // P5 = (A11 + A12) B22
    P = strassen(A11 + A12, B22);
    // subtract P5 from C11 and add it to C21
    C.block_subtract({0, 0}, P);
    C.block_add({0, new_size}, P);

    // P6 = (-A11 + A21)(B11 + B12)
    P = strassen(A21 - A11, B11 + B12);
    // add P6 to C22
    C.block_add({new_size, new_size}, P);

    // P7 = (A12 - A22)(B21 + B22)
    P = strassen(A12 - A22, B21 + B22);
    // add P7 to C11
    C.block_add({0, 0}, P);

    return C;
}

// strassen using dynamic peeling
template<class Scalar>
Matrix<Scalar> multiply_strassen_dynamic(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    // dimension check
    assert(A.cols == B.rows);

    // if any of the dimensions is too small, strassen's algorithm wont help
    if (std::min(A.rows, std::min(A.cols, B.cols)) <= strassen_threshold) {
        return multiply_classic(A, B);
    }

    // split into subblocks:
    // if dimensions are even, make normal strassen step
    // if dimensions are odd, round matrix size down, do a normal strassen algorithm for smaller matrix and perform
    // dynamic peeling after this step finishes
    std::pair<unsigned int, unsigned int>
            block_A = {A.rows / 2, A.cols / 2},
            block_B = {B.rows / 2, B.cols / 2},
            product_block = {A.rows / 2, B.cols / 2};

    Matrix<Scalar> A11 = A.subblock({0, 0}, block_A),
            A12 = A.subblock({0, A.cols / 2}, block_A),
            A21 = A.subblock({A.rows / 2, 0}, block_A),
            A22 = A.subblock({A.rows / 2, A.cols / 2}, block_A);

    Matrix<Scalar> B11 = B.subblock({0, 0}, block_B),
            B12 = B.subblock({0, B.cols / 2}, block_B),
            B21 = B.subblock({B.rows / 2, 0}, block_B),
            B22 = B.subblock({B.rows / 2, B.cols / 2}, block_B);

    // create larger matrix for result
    Matrix<Scalar> C = Matrix<Scalar>::zeros(A.rows, B.cols);

    // temporary matrix, here we will store products
    Matrix<Scalar> P;

    // P1 = (A11 + A22) (B11 + B22)
    P = multiply_strassen_dynamic(A11 + A22, B11 + B22);
    // add P1 to C11 and C22
    C.block_add({0, 0}, P);
    C.block_add(product_block, P);

    // P2 = (A21+ A22) B11
    P = multiply_strassen_dynamic(A21 + A22, B11);
    // add P2 to C21 and subtract from C22
    C.block_add({product_block.first, 0}, P);
    C.block_subtract(product_block, P);

    // P3 = A11 (B12 - B22)
    P = multiply_strassen_dynamic(A11, B12 - B22);
    // add P3 to C12 and C22
    C.block_add({0, product_block.second}, P);
    C.block_add(product_block, P);

    // P4 = A22 (-B11 + B21)
    P = multiply_strassen_dynamic(A22, B21 - B11);
    // add P4 to C11 and C21
    C.block_add({0, 0}, P);
    C.block_add({product_block.first, 0}, P);

    // P5 = (A11 + A12) B22
    P = multiply_strassen_dynamic(A11 + A12, B22);
    // subtract P5 from C11 and add it to C21
    C.block_subtract({0, 0}, P);
    C.block_add({0, product_block.second}, P);

    // P6 = (-A11 + A21)(B11 + B12)
    P = multiply_strassen_dynamic(A21 - A11, B11 + B12);
    // add P6 to C22
    C.block_add(product_block, P);

    // P7 = (A12 - A22)(B21 + B22)
    P = multiply_strassen_dynamic(A12 - A22, B21 + B22);
    // add P7 to C11
    C.block_add({0, 0}, P);

    // fix remaining row and column if dimensions are odd
    dynamic_peeling(A, B, C, 2, 2, 2);

    return C;

}


#endif //FAST_MATRIX_MULTIPLICATION_MULTIPLY_STRASSEN_HPP
