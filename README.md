# metricsC
The metricsC project contains portable routines purely coded in C for estimating econometric models. The project began as routines coded for the .C call in R, which explains why all the functions are void() and function arguments are pointers. The main goal of this project is to provide low-level implementations of econometric models for intuition, higher performance, and insight into the numerical and computational aspects to econometric and statistical model estimation.

metricsC is a project that includes building blocks and functions for estimating econometric models in C. Most of the procedures are void functions callable in software including R in order to speed up computations. To a large extent possible, we provide test (main) programmes for complete C implementation.

The folders organise the programmes into categories/models. The _test.c files contain programmes for implementing the models in C. Also, each  _test.c file contains gcc terminal commands (as comments) for compiling and executing the programmes. For different compilers, a little adaption of the commands may be required.


