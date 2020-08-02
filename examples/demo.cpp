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
    std::cout << "DEMO: Calculating product of matrices A and B." << std::endl;

    int n = 200, m = 190, k = 210;

    std::cout << "A: random " << n << "x" << k << " matrix" << std::endl;
    std::cout << "B: random " << k << "x" << m << " matrix" << std::endl;

    Matrix<double> A = random_float_matrix(n, m);
    Matrix<double> B = random_float_matrix(m, k);
    Matrix<double> classic, result;

    // set output precision
    std::cout << std::fixed << std::setprecision(4);

    Timer t;
    double time;

    // classic multiplication
    t.start();
    classic = multiply_classic(A, B);
    time = t.time_elapsed();

    std::cout << "Classic multiplication took:                \t" << time << "s." << std::endl;

    // strassen static
    t.start();
    result = multiply_strassen_static(A, B);
    time = t.time_elapsed();
    assert(result == classic);  // check if works

    std::cout << "Strassen multiplication <2, 2, 2> (static): \t" << time << "s" << std::endl;

    // strassen dynamic
    t.start();
    result = multiply_strassen_dynamic(A, B);
    time = t.time_elapsed();
    assert(result == classic);  // check if works

    std::cout << "Strassen multiplication <2, 2, 2> (dynamic):\t" << time << "s" << std::endl;

    // laderman (dynamic)
    t.start();
    result = multiply_laderman(A, B);
    time = t.time_elapsed();
    assert(result == classic);  // check if works

    std::cout << "Laderman multiplication <3, 3, 3>:          \t" << time << "s" << std::endl;

    // bini exact multiplication (slow)
    t.start();
    result = multiply_bini_exact(A, B);
    time = t.time_elapsed();
    assert(result == classic);  // check if works

    std::cout << "Bini multiplication <2, 2, 3> (exact):      \t" << time << "s" << std::endl;

    // bini approximative multiplication
    // (useful only for smaller matrices)
    t.start();
    result = multiply_bini(A, B, 1e-6);
    time = t.time_elapsed();

    std::cout << "Bini multiplication <2, 2, 3> (approx):      \t" << time << "s" << std::endl;

    // schonhage exact multiplication (slow)
    t.start();
    result = multiply_schonhage_exact(A, B);
    time = t.time_elapsed();
    assert(result == classic);  // check if works

    std::cout << "Schonhage multiplication <3, 3, 3> (exact):   \t" << time << "s" << std::endl;

    // schonhage approximative multiplication
    // (useful only for smaller matrices)
    t.start();
    result = multiply_bini(A, B, 1e-4);
    time = t.time_elapsed();

    std::cout << "Schonhage multiplication <3, 3, 3> (approx):   \t" << time << "s" << std::endl;

    return 0;
}
