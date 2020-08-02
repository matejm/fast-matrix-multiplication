#ifndef FAST_MATRIX_MULTIPLICATION_MULTIPLY_SCHONHAGE_HPP
#define FAST_MATRIX_MULTIPLICATION_MULTIPLY_SCHONHAGE_HPP

#include <cassert>
#include <vector>
#include "matrix.hpp"
#include "dynamic_peeling.hpp"
#include "polynomial.hpp"

const unsigned int schonhage_threshold = 200;

// Schonhage: Partial and total matrix multiplication,
// Example 2.1
// border rank (<3, 3, 3>) <= 21

// Exact multiplication, problem is that multiplication of polynomials is not O(1) anymore.
template<typename Scalar>
Matrix<Scalar> multiply_schonhage_exact(const Matrix<Scalar> &A, const Matrix<Scalar> &B) {
    // convert to polynomials
    Matrix<Polynomial<Scalar>> poly_A(A);
    Matrix<Polynomial<Scalar>> poly_B(B);

    Matrix<Polynomial<Scalar>> poly_C = multiply_schonhage(poly_A, poly_B, Polynomial<Scalar>::epsilon());

    // convert back from polynomials
    return polynomial_to_scalar(poly_C);
}

// actual Schonhage's algorithm
template<typename Poly>
Matrix<Poly> multiply_schonhage(const Matrix<Poly> &_A, const Matrix<Poly> &_B, const Poly &epsilon) {
    // check dimensions
    assert(_A.cols == _B.rows);

    // if matrices are too small for Bini's algorithm we have nothing to do but multiply it classicaly
    if (_A.rows < 3 || _A.cols < 3 || _B.cols < 3) {
        return multiply_classic(_A, _B);
    } else if (_A.rows <= schonhage_threshold || _A.cols <= schonhage_threshold || _B.cols <= schonhage_threshold) {
        return multiply_classic(_A, _B);
    }

    // create subblocks
    // last line and up to 2 last columns may not be included, this is handled by dynamic peeling
    unsigned int block_rows_A = _A.rows / 3,
            block_cols_A = _A.cols / 3,
            block_rows_B = _B.rows / 3,
            block_cols_B = _B.cols / 3;

    std::pair<unsigned int, unsigned int>
            block_A = {block_rows_A, block_cols_A},
            block_B = {block_rows_B, block_cols_B};

    // | A11 A12 A13 | | B11 B12 B13 |
    // | A21 A22 A23 | | B21 B22 B23 |
    // | A31 A32 A33 | | B31 B32 B33 |

    // matrices in blocks
    std::vector<std::vector<Matrix<Poly>>>
            A(3, std::vector<Matrix<Poly>>(3)),
            B(3, std::vector<Matrix<Poly>>(3)),
            C(3, std::vector<Matrix<Poly>>(3));

    // subtract 1 because indices start at 1
    A[1 - 1][1 - 1] = _A.subblock({0, 0}, block_A),
    A[1 - 1][2 - 1] = _A.subblock({0, block_cols_A}, block_A),
    A[1 - 1][3 - 1] = _A.subblock({0, 2 * block_cols_A}, block_A),
    A[2 - 1][1 - 1] = _A.subblock({block_rows_A, 0}, block_A),
    A[2 - 1][2 - 1] = _A.subblock({block_rows_A, block_cols_A}, block_A),
    A[2 - 1][3 - 1] = _A.subblock({block_rows_A, 2 * block_cols_A}, block_A),
    A[3 - 1][1 - 1] = _A.subblock({2 * block_rows_A, 0}, block_A),
    A[3 - 1][2 - 1] = _A.subblock({2 * block_rows_A, block_cols_A}, block_A),
    A[3 - 1][3 - 1] = _A.subblock({2 * block_rows_A, 2 * block_cols_A}, block_A);

    B[1 - 1][1 - 1] = _B.subblock({0, 0}, block_B),
    B[1 - 1][2 - 1] = _B.subblock({0, block_cols_B}, block_B),
    B[1 - 1][3 - 1] = _B.subblock({0, 2 * block_cols_B}, block_B),
    B[2 - 1][1 - 1] = _B.subblock({block_rows_B, 0}, block_B),
    B[2 - 1][2 - 1] = _B.subblock({block_rows_B, block_cols_B}, block_B),
    B[2 - 1][3 - 1] = _B.subblock({block_rows_B, 2 * block_cols_B}, block_B),
    B[3 - 1][1 - 1] = _B.subblock({2 * block_rows_B, 0}, block_B),
    B[3 - 1][2 - 1] = _B.subblock({2 * block_rows_B, block_cols_B}, block_B),
    B[3 - 1][3 - 1] = _B.subblock({2 * block_rows_B, 2 * block_cols_B}, block_B);

    // create new empty matrix for product
    // | C11 C12 C13 |
    // | C21 C22 C23 |
    // | C31 C32 C33 |
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            // fill with zeros
            C[i][j] = Matrix<Poly>::zeros(block_rows_A, block_cols_B);
        }
    }

    // epsilon squared
    Poly epsilon2 = epsilon * epsilon;

    // actual algorithm,
    // see article Partial and Total Matrix Multiplication, example 2.2 for formulas
    for (unsigned int i = 0; i < 3; i++) {
        // Wi
        Matrix<Poly> W = multiply_schonhage(
                A[i][0],
                B[1][i] + B[2][i],
                epsilon
        );

        // D'ji from article, we are calculating C
        // Cji = 1/epsilon^2 (Uij + Vij - Wi) + 1/epsilon (Vji - Vjj)
        // Cii = 1/epsilon^2 (Uii + Vii - Wi)
        // because we do not want to keep all matrices in memory,
        // this calculation will be split in multiple parts
        for (unsigned int j = 0; j < 3; j++) {
            Matrix<Poly> U, V;

            if (j == i) {
                // Uii
                U = multiply_schonhage(
                        A[i][0] + epsilon2 * A[i][1],
                        epsilon2 * B[0][i] + B[1][i],
                        epsilon
                );
                // Vii
                V = multiply_schonhage(
                        A[i][0] + epsilon2 * A[i][2],
                        B[2][i],
                        epsilon
                );

                Matrix<Poly> to_subtract = V;
                to_subtract /= epsilon;
                // subtract Vjj from all Cji (where i != j)
                for (unsigned int k = 0; k < 3; k++) {
                    if (k != j) {
                        C[j][k] -= to_subtract;
                    }
                }

            } else {
                // Uij
                U = multiply_schonhage(
                        A[i][0] + epsilon2 * A[j][1],
                        B[1][i] - epsilon * B[0][j],
                        epsilon
                );
                // Vij
                V = multiply_schonhage(
                        A[i][0] + epsilon2 * A[j][2],
                        B[2][i] + epsilon * B[0][j],
                        epsilon
                );
                // 1/epsilon Vji
                Matrix<Poly> to_add = V;
                to_add /= epsilon;
                C[i][j] += to_add;
            }

            // 1/epsilon^2 (Uij + Vij - Wi)
            // (or if i == j (Uii + Vii - Wi)
            Matrix<Poly> to_add = U + V - W;
            to_add /= epsilon2;
            C[j][i] += to_add;
        }
    }

    // join block from matrix C into a single block
    Matrix<Poly> _C = Matrix<Poly>::zeros(_A.rows, _B.cols);

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            _C.block_add({block_rows_A * i, block_cols_B * j}, C[i][j]);
        }
    }

    // dynamic peeling for not included rows and cols
    dynamic_peeling(_A, _B, _C, 3, 3, 3);

    return _C;
}


#endif //FAST_MATRIX_MULTIPLICATION_MULTIPLY_SCHONHAGE_HPP
