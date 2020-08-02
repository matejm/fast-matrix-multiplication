#ifndef FAST_MATRIX_MULTIPLICATION_POLYNOMIAL_HPP
#define FAST_MATRIX_MULTIPLICATION_POLYNOMIAL_HPP

#include "matrix.hpp"
#include <iostream>
#include <vector>

// Represents polynomial in epsilon, needed for Bini's algorithm.
// Because algorithm requires polynomial multiplication, which is quadratic, this is not a great idea.

template<typename Scalar>
class Polynomial {
public:
    std::vector<Scalar> a;

    // linear polynomial
    Polynomial(Scalar a0, Scalar a1) : a({a0, a1}) {};

    // Returns a constant polynomial with a constant c.
    Polynomial(Scalar a0) : a({a0}) {};

    // zero
    Polynomial() : a({Scalar(0)}) {}

    // Returns a polynomial that represents epsilon
    static Polynomial<Scalar> epsilon() {
        return Polynomial<Scalar>(Scalar(0), Scalar(1));
    }

    Polynomial<Scalar> &operator*=(const Polynomial<Scalar> &other) {
        std::vector<Scalar> product(a.size() + other.a.size() - 1, 0);

        for (size_t i = 0; i < other.a.size(); ++i) {
            for (size_t j = 0; j < a.size(); ++j) {
                product[i + j] += a[j] * other.a[i];
            }
        }

        a = product;
        return *this;
    }

    Polynomial<Scalar> &operator+=(const Polynomial<Scalar> &other) {
        while (a.size() < other.a.size()) {
            a.push_back(Scalar(0));
        }
        for (unsigned int i = 0; i < std::min(a.size(), other.a.size()); i++) {
            a[i] += other.a[i];
        }
        return *this;
    }

    Polynomial<Scalar> &operator-=(const Polynomial<Scalar> &other) {
        while (a.size() < other.a.size()) {
            a.push_back(Scalar(0));
        }
        for (unsigned int i = 0; i < std::min(a.size(), other.a.size()); i++) {
            a[i] -= other.a[i];
        }
        return *this;
    }

    Polynomial<Scalar> &operator/=(const Polynomial<Scalar> &other) {
        // Not actual polynomial division.
        // Needed for Bini's algorithm where we divide with epsilon at the end of a recursive step.
        // Needed for Schonhage's algorithm where we divide with epsilon and epsilon squared.
        assert(// epsilon
                (other.a.size() == 2 && other.a[0] == Scalar(0) && other.a[1] == Scalar(1)) ||
                // epsilon squared
                (other.a.size() == 3 && other.a[0] == Scalar(0) && other.a[1] == Scalar(0) && other.a[2] == Scalar(1))
        );

        // divide by epsilon
        if (other.a.size() == 2) {
            if (a.size() == 1) {
                a = {0};
            } else {
                for (unsigned int i = 0; i < a.size() - 1; i++) {
                    a[i] = a[i + 1];
                }
                a.pop_back();
            }
        } else if (other.a.size() == 3) {
            if (a.size() <= 2) {
                a = {0};
            } else {
                for (unsigned int i = 0; i < a.size() - 2; i++) {
                    a[i] = a[i + 2];
                }
                a.pop_back();
                a.pop_back();
            }
        } else {
            assert("Not implemented.");
        }

        return *this;
    }
};

template<class Scalar>
Polynomial<Scalar> operator+(const Polynomial<Scalar> &p, const Polynomial<Scalar> &q) {
    Polynomial<Scalar> result = p;
    result += q;
    return result;
}

template<class Scalar>
Polynomial<Scalar> operator-(const Polynomial<Scalar> &p, const Polynomial<Scalar> &q) {
    Polynomial<Scalar> result = p;
    result -= q;
    return result;
}

template<class Scalar>
Polynomial<Scalar> operator*(const Polynomial<Scalar> &p, const Polynomial<Scalar> &q) {
    Polynomial<Scalar> result = p;
    result *= q;
    return result;
}

template<class Scalar>
Polynomial<Scalar> operator/(const Polynomial<Scalar> &p, const Polynomial<Scalar> &q) {
    Polynomial<Scalar> result = p;
    result /= q;
    return result;
}

// function that converts polynomial matrix to scalar matrix
// (we can imagine that we sent epsilon to zero, actually we just take constant term)
// to convert in other direction just use matrix constructor
template<typename Scalar>
Matrix<Scalar> polynomial_to_scalar(const Matrix<Polynomial<Scalar>> &A) {
    std::vector<Scalar> data(A.data.size());

    for (unsigned int i = 0; i < A.data.size(); i++) {
        data[i] = (A.data[i]).a[0];
    }

    return Matrix<Scalar>(data, A.rows, A.cols);
}

// make poylnomial printable
template<typename Scalar>
std::ostream &operator<<(std::ostream &os, const Polynomial<Scalar> &p) {
    for (unsigned int i = 0; i < p.a.size() - 1; i++) {
        os << p.a[i] << "e^" << i << " + ";
    }
    return os << p.a[p.a.size() - 1] << "e^" << p.a.size() - 1;
}

#endif //FAST_MATRIX_MULTIPLICATION_POLYNOMIAL_HPP
