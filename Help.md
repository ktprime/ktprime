Ktprime Project Introduction
Overview
The ktprime project is a collection of useful programs, including prime/twin sieves, hash maps, timers, LRU caches, and top - k algorithms. This repository aims to provide efficient and practical algorithms and data structures for various computational tasks, especially those related to prime number calculations.
Key Components and Their Functions
1. Prime Number Calculation
PrimeNumber.cpp:
allocWheelBlock: Allocates memory for wheel blocks used in prime number sieve algorithms. It initializes the wheel pool and stock cache for efficient prime number sieving.
initMediumWheel: Initializes the medium wheel sieve, which is used to pre - process prime numbers in a specific range. It marks non - prime numbers based on the wheel sieve principle.
Ktprime.cpp:
initPattern: Initializes the prime number patterns. It calculates the number of patterns and fills the pattern array according to different configuration parameters.
initMoudle: Initializes the module for prime number calculations. It calculates the inverse of the wheel modulo each prime number and stores the results in the module array.
executeCmd: Parses and executes commands related to prime number calculations, such as benchmarking, unit testing, and calculating prime k - tuples.
TwinPrime.cpp:
crossMedium1 and crossMedium3: These functions are used to cross out non - prime numbers in the medium sieve algorithm. They mark multiples of prime numbers in the bit array.
eratSieveMedium: Performs the medium sieve algorithm to find prime numbers in a given range.
executeCmd: Parses and executes commands for twin prime calculations, including benchmarking and sieving in specific ranges.
SophieGermain.cpp:
getPattern: Generates prime number patterns based on the wheel for Sophie Germain prime calculations.
executeCmd: Parses and executes commands for Sophie Germain prime calculations, such as benchmarking, listing differences, and calculating primes in specific ranges.
2. Hash Table and LRU Cache
ehash_table.hpp:
CTZ: Counts the leading zero bits of a number, which is used in hash table operations.
reserve: Reserves space for a given number of elements in the hash table. If necessary, it rehashes the table to accommodate more elements.
operator[]: Inserts or retrieves an element from the hash table. If the key does not exist, it inserts a new element.
ehash_lru.h:
nextid: Generates a unique identifier for each entry in the LRU cache.
reserve: Reserves space for a given number of elements in the LRU cache. It may remove half of the elements if the cache is full.
operator[]: Inserts or retrieves an element from the LRU cache. If the key does not exist, it inserts a new element and updates the LRU order.
3. Other Utilities
abx.cpp: Implements a calculator for linear equations. It provides a Fraction class to handle fractional arithmetic and a Calculator class to calculate the result of input expressions.
wheel_timer.h: Defines an interface TaskRunnerInterface for task scheduling and cancellation. It provides functions to schedule tasks to run after a given delay and cancel scheduled tasks.
FastGn.cpp: Provides a fast algorithm for calculating Goldbach partitions. It uses segmented sieve methods to improve performance. The executeCmd function parses and executes commands related to Goldbach partition calculations.
How to Use
To use these programs, you can include the relevant header files in your project and call the corresponding functions. For example, if you want to use the hash table, you can include ehash_table.hpp and use the provided functions such as operator[] and try_get.

For prime number calculations, you can call the executeCmd functions in Ktprime.cpp, TwinPrime.cpp, etc., with appropriate commands to perform different operations such as benchmarking, unit testing, and calculating prime numbers in specific ranges.
License
Please refer to the project's license file for detailed licensing information.
Contribution
Contributions to this project are welcome. If you have any improvements or new features to add, please submit a pull request.
