# CMC MSU course "Distributed systems"

The course contained two tasks:

1. Develop a program that implements the specified algorithm and obtain a time estimation of its performance.

2. Refine the MPI program implemented as part of the "Supercomputers and Parallel Data Processing" course. Using parallel I/O (MPI-IO), add checkpoints to allow the program to continue in the event of a failure. Prepare a report on the completed assignment, including a description of the algorithm, implementation details, and time performance evaluations. The report and the program code must be posted in the tracker. Implement one of the following three scenarios for post-failure operation:
    a) continue execution only on the “healthy” (non-failed) processes;
    b) instead of the failed processes, create new MPI processes to continue the calculations;
    c) when starting the program, immediately launch an additional number of MPI processes that can be used in case of a failure.
