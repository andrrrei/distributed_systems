# CMC MSU Course "Distributed Systems"

The course included two tasks:

## Task 1: Algorithm Implementation and Performance Evaluation

Develop a program that implements the specified algorithm and obtain a time estimation of its performance.

## Task 2: Enhancing an MPI Program with Fault Tolerance

Refine the MPI program implemented as part of the **"Supercomputers and Parallel Data Processing"** course. Using **parallel I/O (MPI-IO)**, add checkpoints to allow the program to continue in the event of a failure. Implement one of the following three scenarios for post-failure operation:

1. **Continue on healthy processes:** Continue execution only on the “healthy” (non-failed) processes.
2. **Replace failed processes:** Create new MPI processes to replace those that failed, and use these new processes to continue the calculations.
3. **Use spare processes:** When starting the program, immediately launch an additional number of MPI processes that can be used in case of a failure.

### Deliverables:
- Prepare a **PDF report** describing the solution, including:
  - Algorithm description.
  - Implementation details.
  - Performance evaluation.
- Submit the report and program source code in the designated folder.

Each task is placed in its corresponding folder and contains a PDF report with a description of the solution.
