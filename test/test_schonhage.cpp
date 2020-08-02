#include <gtest/gtest.h>

#include "helpers.hpp"

#include "multiply_schonhage.hpp"
#include "multiply_classic.hpp"
#include "matrix.hpp"

TEST(SchonhageExact, Basic) {
    // few basic tests that can be calculated by hand
    Matrix<int> A, B;

    // too small for algorithm
    A = Matrix<int>({1, 2, 3, 4}, 2, 2);
    B = Matrix<int>({4, 3, 2, 1}, 2, 2);

    ASSERT_EQ(
            multiply_schonhage_exact(A, B), Matrix<int>({8, 5, 20, 13}, 2, 2)
    );

    // scalar
    A = Matrix<int>(1, 1, 2);
    B = Matrix<int>(1, 1, 3);

    ASSERT_EQ(
            multiply_schonhage_exact(A, B), Matrix<int>(1, 1, 6)
    );

    // basic algorithm step
    A = Matrix<int>({1, 2, 3, 4, 5, 6, 7, 8, 9}, 3, 3);
    B = Matrix<int>({9, 8, 7, 6, 5, 4, 3, 2, 1}, 3, 3);

    ASSERT_EQ(
            multiply_schonhage_exact(A, B), Matrix<int>({30, 24, 18, 84, 69, 54, 138, 114, 90}, 3, 3)
    );

    // basic step and dynamic peeling
    A = Matrix<int>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, 3, 4);
    B = Matrix<int>({12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1}, 4, 3);

    ASSERT_EQ(
            multiply_schonhage_exact(A, B), Matrix<int>({60, 50, 40, 180, 154, 128, 300, 258, 216}, 3, 3)
    );
}

TEST(SchonhageExact, Square) {
    // Test random square matrices from size 1x1 to 50x50
    Matrix<int> A, B;

    for (unsigned int i = 1; i <= 50; i++) {
        A = random_int_matrix(i, i);
        B = random_int_matrix(i, i);

        ASSERT_EQ(
                multiply_classic(A, B),
                multiply_schonhage_exact(A, B)
        );
    }
}

TEST(SchonhageExact, NonSquare) {
    // Test 50 random matrices of 'random' size,
    // size is at max 50.
    unsigned int n, k, m;
    Matrix<int> A, B;

    for (unsigned int i = 0; i < 50; i++) {
        n = rand() % 50 + 1;
        k = rand() % 50 + 1;
        m = rand() % 50 + 1;

        A = random_int_matrix(n, k);
        B = random_int_matrix(k, m);

        ASSERT_EQ(
                multiply_classic(A, B),
                multiply_schonhage_exact(A, B)
        );
    }
}

TEST(SchonhageExact, Large) {
    // Test a few larger matrices
    Matrix<int> A, B;

    A = random_int_matrix(111, 111);
    B = random_int_matrix(111, 111);
    ASSERT_EQ(multiply_classic(A, B), multiply_schonhage_exact(A, B));

    A = random_int_matrix(100, 123);
    B = random_int_matrix(123, 100);
    ASSERT_EQ(multiply_classic(A, B), multiply_schonhage_exact(A, B));

    A = random_int_matrix(123, 321);
    B = random_int_matrix(321, 21);
    ASSERT_EQ(multiply_classic(A, B), multiply_schonhage_exact(A, B));
}


TEST(SchonhageApprox, Basic) {
    // few basic tests that can be calculated by hand
    Matrix<double> A, B;

    // too small for algorithm
    A = Matrix<double>({1, 2, 3, 4}, 2, 2);
    B = Matrix<double>({4, 3, 2, 1}, 2, 2);

    ASSERT_LE(
            maximum_relative_difference(
                    multiply_classic(A, B),
                    multiply_schonhage(A, B, 1e-6)),
            0.01
    );

    // scalar
    A = Matrix<int>(1, 1, 2);
    B = Matrix<int>(1, 1, 3);

    ASSERT_LE(
            maximum_relative_difference(
                    multiply_classic(A, B),
                    multiply_schonhage(A, B, 1e-6)),
            0.01
    );

    // basic algorithm step
    A = Matrix<int>({1, 2, 3, 4, 5, 6, 7, 8, 9}, 3, 3);
    B = Matrix<int>({9, 8, 7, 6, 5, 4, 3, 2, 1}, 3, 3);

    ASSERT_LE(
            maximum_relative_difference(
                    multiply_classic(A, B),
                    multiply_schonhage(A, B, 1e-6)),
            0.01
    );

    // basic step and dynamic peeling
    A = Matrix<int>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, 3, 4);
    B = Matrix<int>({12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1}, 4, 3);

    ASSERT_LE(
            maximum_relative_difference(
                    multiply_classic(A, B),
                    multiply_schonhage(A, B, 1e-6)),
            0.01
    );
}

TEST(SchonhageApprox, Square) {
    // Test random square matrices from size 1x1 to 50x50
    Matrix<double> A, B;

    for (unsigned int i = 1; i <= 50; i++) {
        A = random_float_matrix(i, i);
        B = random_float_matrix(i, i);

        ASSERT_LE(
                maximum_relative_difference(
                        multiply_classic(A, B),
                        multiply_schonhage(A, B, 1e-2)),
                0.1
        );
    }
}

TEST(SchonhageApprox, NonSquare) {
    // Test 50 random matrices of 'random' size,
    // size is at max 50.
    unsigned int n, k, m;
    Matrix<double> A, B;

    for (unsigned int i = 0; i < 50; i++) {
        n = rand() % 50 + 1;
        k = rand() % 50 + 1;
        m = rand() % 50 + 1;

        A = random_float_matrix(n, k);
        B = random_float_matrix(k, m);

        ASSERT_LE(
                maximum_relative_difference(
                        multiply_classic(A, B),
                        multiply_schonhage(A, B, 1e-2)),
                0.1
        );
    }
}

TEST(SchonhageApprox, Large) {
    // Test a few larger matrices
    Matrix<double> A, B;

    A = random_float_matrix(111, 111);
    B = random_float_matrix(111, 111);

    ASSERT_LE(
            maximum_relative_difference(
                    multiply_classic(A, B),
                    multiply_schonhage(A, B, 1e-1)),
            0.1
    );

    A = random_float_matrix(100, 123);
    B = random_float_matrix(123, 100);

    ASSERT_LE(
            maximum_relative_difference(
                    multiply_classic(A, B),
                    multiply_schonhage(A, B, 1e-1)),
            0.1
    );

    A = random_float_matrix(123, 321);
    B = random_float_matrix(321, 21);

    ASSERT_LE(
            maximum_relative_difference(
                    multiply_classic(A, B),
                    multiply_schonhage(A, B, 1e-1)),
            0.1
    );
}