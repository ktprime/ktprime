README
Project Overview
This repository, ktprime, contains a collection of useful programs, including prime and twin sieve algorithms, hash map implementations, timer utilities, LRU caches, and top-k algorithms. The following is a detailed introduction to the main functions of each file:
File Descriptions
1. Ktprime.cpp
Functionality:
Prime Counting: It provides functions to count prime numbers, twin primes, and ktuplet primes. The main counting functions are Ktprime which can calculate the number of primes based on different patterns.
Command Parsing: The parseCmd function parses command-line parameters and sets corresponding configuration flags according to the input commands.
Initialization and Main Execution: The initKtprime function initializes the program, and the main function is the entry point, which handles command-line arguments and user input in a loop.
Printing Information: Functions like printInfo and printProgress are used to print program information and progress during execution.
2. TwinPrime.cpp
Functionality:
Fast Segmented Sieving: It implements a fast segmented sieving algorithm for twin primes, which is considered one of the fastest algorithms before 2019.
Benchmarking: It provides benchmark data for different ranges of numbers, showing the performance of the algorithm on various CPU models.
Command Parsing: The parseCmd function parses command-line parameters and performs corresponding operations such as setting sieve size, cache size, and cache segments.
3. topK.cpp
Functionality:
Top-K Algorithms: It contains multiple algorithms to find the top-k elements in an array, including STL sorting, nth_element, priority queue, k-ary heap, min-max heap, and d-ary heap algorithms.
Performance Testing: The program measures the execution time of each algorithm and checks the results to ensure correctness.
4. LICENSE
Functionality:
It is an MIT license file, which grants users the right to use, copy, modify, merge, publish, distribute, sublicense, and/or sell the software, subject to the conditions specified in the license.
Usage
Compilation: Compile the source code using a C++ compiler with appropriate flags, such as -DSIEVE_SIZE=2048 -march=native -funroll-loops -O3 -s -pipe.
Execution: Run the compiled executable and pass command-line parameters according to the requirements. For example, you can specify the range of numbers to count primes or set configuration options.
Dependencies
Some of the code may rely on specific libraries, such as pdqsort.h, minmax_and_dary_heap.hpp, and ska_sort.hpp. Make sure these libraries are available in the include path if you want to use the corresponding features.
Note
The code is mainly written in C++, and some parts may use platform-specific features. Please adjust the compilation and execution environment according to your system.
Contributing
Contributions to this project are welcome. If you have any improvements or bug fixes, please submit a pull request.
