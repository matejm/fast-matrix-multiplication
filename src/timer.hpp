#ifndef FAST_MATRIX_MULTIPLICATION_TIMER_HPP
#define FAST_MATRIX_MULTIPLICATION_TIMER_HPP

#include <chrono>

// simple class for measuring elapsed time
class Timer {
public:
    typedef std::chrono::high_resolution_clock Clock;

    void start() {
        epoch = Clock::now();
    }

    // time elapsed in nanoseconds
    Clock::duration time_elapsed_exact() const {
        return Clock::now() - epoch;
    }

    // time elapsed in seconds
    double time_elapsed() const {
        return std::chrono::duration_cast<std::chrono::microseconds>(time_elapsed_exact()).count() / 1000000.0;
    }


private:
    Clock::time_point epoch;
};

#endif //FAST_MATRIX_MULTIPLICATION_TIMER_HPP
