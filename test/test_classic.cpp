#include <gtest/gtest.h>

#include "helpers.hpp"

#include "multiply_classic.hpp"
#include "matrix.hpp"

TEST(Classic, Basic) {
    // few basic tests that can be calculated by hand
    Matrix<int> A, B;

    // <2, 2, 2>
    A = Matrix<int>({1, 2, 3, 4}, 2, 2);
    B = Matrix<int>({4, 3, 2, 1}, 2, 2);

    ASSERT_EQ(
            multiply_classic(A, B),
            Matrix<int>({8, 5, 20, 13}, 2, 2)
    );

    // <2, 2, 3>
    A = Matrix<int>({1, 2, 3, 4}, 2, 2);
    B = Matrix<int>({6, 5, 4, 3, 2, 1}, 2, 3);

    ASSERT_EQ(
            multiply_classic(A, B),
            Matrix<int>({12, 9, 6, 30, 23, 16}, 2, 3)
    );

    // scalar multiplication test
    A = Matrix<int>(1, 1, 2);
    B = Matrix<int>(1, 1, 3);

    ASSERT_EQ(
            multiply_classic(A, B),
            Matrix<int>(1, 1, 6)
    );

    A = Matrix<int>({1, 2, 3, 4, 5, 6, 7, 8}, 2, 4);
    B = Matrix<int>({8, 7, 6, 5, 4, 3, 2, 1}, 4, 2);

    ASSERT_EQ(
            multiply_classic(A, B),
            Matrix<int>({40, 30, 120, 94}, 2, 2)
    );

    // a bit larger matrix
    A = Matrix<int>({1, 2, 3, 4, 5, 6, 7, 8}, 2, 4);
    B = Matrix<int>({12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1}, 4, 3);

    ASSERT_EQ(
            multiply_classic(A, B),
            Matrix<int>({60, 50, 40, 180, 154, 128}, 2, 3)
    );
}

TEST(Classic, Idenity) {
    // really trivial test, but still better than nothing
    Matrix<int> A, A_squared, A_cubed;

    for (unsigned int size = 1; size < 50; size++) {
        std::vector<int> data(size * size);

        for (unsigned int i = 0; i < size; ++i) {
            data[i * size + i] = 1;
        }

        A = Matrix<int>(data, size, size);

        A_squared = multiply_classic(A, A);
        A_cubed = multiply_classic(A_squared, A);

        ASSERT_EQ(A, A_squared);
        ASSERT_EQ(A, A_cubed);
    }
}