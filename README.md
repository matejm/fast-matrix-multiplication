## Hitro množenje matrik

Ta repozitorij vsebuje implementacijo izbranih algoritmov množenja matrik in je del moje diplomske naloge. 

## Fast Matrix Multiplication

This repository contains implementation of next fast matrix multiplication algorithms:
- classic matrix multiplication
- Strassen's algorithm (`<2, 2, 2>`, tensor rank 7)
- Laderman's algorithm (`<2, 2, 2>`, tensor rank 23)
- Bini's approximative algorithm (`<2, 2, 3>`, tensor rank 10)
- Schonhage's approximative algorithm (`<3, 3, 3>`, tensor rank 21)

Algorithms are far from optimized, the main goal was to see the diffenences in time complexitiy.

#### Instalation and running examples

Clone this repository
```
git clone git@github.com:matejm/fast-matrix-multiplication.git
```
Go to root folder of the repository and do the usual compile & run cycle:
```
mkdir -p build
cd build
cmake ..
make
```
Run demo of all algorithms or time algorithm execution for different matrix sizes:
```
./../bin/demo
./../bin/time_multiplication
```

#### Testing

Google's testing framework is used. To install it on Ubuntu, run:
```
sudo apt install libgtest-dev
```
To test the algorithms, run:
```
make tests
./../bin/tests
```

