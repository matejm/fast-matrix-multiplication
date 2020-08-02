#include <gtest/gtest.h>

#include "helpers.hpp"

#include "multiply_strassen.hpp"
#include "multiply_classic.hpp"
#include "matrix.hpp"

TEST(StrassenStatic, Basic) {
    // few basic tests that can be calculated by hand
    Matrix<int> A, B;

    A = Matrix<int>({1, 2, 3, 4}, 2, 2);
    B = Matrix<int>({4, 3, 2, 1}, 2, 2);

    ASSERT_EQ(
            multiply_strassen_static(A, B), Matrix<int>({8, 5, 20, 13}, 2, 2)
    );

    A = Matrix<int>(1, 1, 2);
    B = Matrix<int>(1, 1, 3);

    ASSERT_EQ(
            multiply_strassen_static(A, B), Matrix<int>(1, 1, 6)
    );

    A = Matrix<int>({1, 2, 3, 4, 5, 6, 7, 8}, 2, 4);
    B = Matrix<int>({8, 7, 6, 5, 4, 3, 2, 1}, 4, 2);

    ASSERT_EQ(
            multiply_strassen_static(A, B), Matrix<int>({40, 30, 120, 94}, 2, 2)
    );
}

TEST(StrassenStatic, Square) {
    // Test random square matrices from size 1x1 to 50x50
    Matrix<int> A, B;

    for (unsigned int i = 1; i <= 50; i++) {
        A = random_int_matrix(i, i);
        B = random_int_matrix(i, i);

        ASSERT_EQ(
                multiply_classic(A, B), multiply_strassen_static(A, B)
        );
    }
}

TEST(StrassenStatic, NonSquare) {
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
                multiply_classic(A, B), multiply_strassen_static(A, B)
        );
    }
}

TEST(StrassenStatic, Large) {
    // Test a few larger matrices
    Matrix<int> A, B;

    A = random_int_matrix(111, 111);
    B = random_int_matrix(111, 111);
    ASSERT_EQ(multiply_classic(A, B), multiply_strassen_static(A, B));

    A = random_int_matrix(100, 123);
    B = random_int_matrix(123, 100);
    ASSERT_EQ(multiply_classic(A, B), multiply_strassen_static(A, B));

    A = random_int_matrix(123, 321);
    B = random_int_matrix(321, 21);
    ASSERT_EQ(multiply_classic(A, B), multiply_strassen_static(A, B));
}


TEST(StrassenDynamic, Basic) {
    // few basic tests that can be calculated by hand
    Matrix<int> A, B;

    A = Matrix<int>({1, 2, 3, 4}, 2, 2);
    B = Matrix<int>({4, 3, 2, 1}, 2, 2);

    ASSERT_EQ(
            multiply_strassen_dynamic(A, B), Matrix<int>({8, 5, 20, 13}, 2, 2)
    );

    A = Matrix<int>(1, 1, 2);
    B = Matrix<int>(1, 1, 3);

    ASSERT_EQ(
            multiply_strassen_dynamic(A, B), Matrix<int>(1, 1, 6)
    );

    A = Matrix<int>({1, 2, 3, 4, 5, 6, 7, 8}, 2, 4);
    B = Matrix<int>({8, 7, 6, 5, 4, 3, 2, 1}, 4, 2);

    ASSERT_EQ(
            multiply_strassen_dynamic(A, B), Matrix<int>({40, 30, 120, 94}, 2, 2)
    );
}

TEST(StrassenDynamic, Square) {
    // Test random square matrices from size 1x1 to 50x50
    Matrix<int> A, B;

    for (unsigned int i = 1; i <= 50; i++) {
        A = random_int_matrix(i, i);
        B = random_int_matrix(i, i);

        ASSERT_EQ(
                multiply_classic(A, B), multiply_strassen_dynamic(A, B)
        );
    }
}

TEST(StrassenDynamic, NonSquare) {
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
                multiply_classic(A, B), multiply_strassen_dynamic(A, B)
        );
    }
}

TEST(StrassenDynamic, Large) {
    // Test a few larger matrices
    Matrix<int> A, B;

    A = random_int_matrix(111, 111);
    B = random_int_matrix(111, 111);
    ASSERT_EQ(multiply_classic(A, B), multiply_strassen_dynamic(A, B));

    A = random_int_matrix(100, 123);
    B = random_int_matrix(123, 100);
    ASSERT_EQ(multiply_classic(A, B), multiply_strassen_dynamic(A, B));

    A = random_int_matrix(123, 321);
    B = random_int_matrix(321, 21);
    ASSERT_EQ(multiply_classic(A, B), multiply_strassen_dynamic(A, B));
}
