#include <iomanip>
#include <iostream>
#include "matrix.hpp"
#include "util.hpp"
#include "multiply_classic.hpp"
#include "multiply_strassen.hpp"
#include "multiply_bini.hpp"
#include "multiply_laderman.hpp"
#include "multiply_schonhage.hpp"
#include "timer.hpp"

int main() {
    // max matrix size
    unsigned int
            N = 3000,
            step = 100,
            start = 100;

    std::cout << "Timing matrix product calculation." << std::endl;
    std::cout << "Each algorithm is tested on matrix sizes from " << start << " to " << N << ", step " << step << "."
              << std::endl;

    // set output precision
    std::cout << std::fixed << std::setprecision(4);

    Timer t;
    double time;
    Matrix<double> classic, result;

    unsigned int size = start;

    while (size <= N) {
        std::cout << std::endl << "SIZE " << size << "x" << size << std::endl;

        // generate random matrices
        // (this test uses only square matrices)
        Matrix<double> A = random_float_matrix(size, size);
        Matrix<double> B = random_float_matrix(size, size);

        // classical multiplication
        t.start();
        classic = multiply_classic(A, B);
        time = t.time_elapsed();

        std::cout << "Classic multiplication took:                \t" << time << "s" << std::endl;

        // strassen static
        t.start();
        result = multiply_strassen_static(A, B);
        time = t.time_elapsed();

       std::cout << "Strassen multiplication <2, 2, 2> (static): \t" << time << "s" << std::endl;

        // strassen dynamic
        t.start();
        result = multiply_strassen_dynamic(A, B);
        time = t.time_elapsed();
        assert(result == classic);

        std::cout << "Strassen multiplication <2, 2, 2> (dynamic):\t" << time << "s" << std::endl;

        // laderman (dynamic)
        t.start();
        result = multiply_laderman(A, B);
        time = t.time_elapsed();
        assert(result == classic);  // check if works

        std::cout << "Laderman multiplication <3, 3, 3>:          \t" << time << "s" << std::endl;

        // bini exact multiplication
        // (slow)
        t.start();
        result = multiply_bini_exact(A, B);
        time = t.time_elapsed();

        std::cout << "Bini multiplication <2, 2, 3> (exact):       \t" << time << "s" << std::endl;

        // bini approximative multiplication
        // (useful only for smaller matrices, for larger matrices it is possible that
        // result is not even close to actual result)
        t.start();
        result = multiply_bini(A, B, 1e-1);
        time = t.time_elapsed();

        std::cout << "Bini multiplication <2, 2, 3> (approx):      \t" << time << "s" << std::endl;

        // bini exact multiplication
        // (slow)
        t.start();
        result = multiply_schonhage_exact(A, B);
        time = t.time_elapsed();

        std::cout << "Schonhage multiplication <3, 3, 3> (exact):    \t" << time << "s" << std::endl;

        // schonhage approximative multiplication
        // (useful only for smaller matrices, for larger matrices it is possible
        // that result is not even close to actual result)
        t.start();
        result = multiply_schonhage(A, B, 1e-1);
        time = t.time_elapsed();

        std::cout << "Schonhage multiplication <3, 3, 3> (approx):   \t" << time << "s" << std::endl;

        size += step;
    }

    return 0;
}

