This project is build and written in Visual Studio Code with CMake extension
To run this project on HPCC, compile a seperate executable using "g++ -Wall -fopenmp -o (executable name) main.cpp"

## Notes about test functions:
> *Test(int rowA, int colA_rowB, int colB, int attempts)
The test functions found in "TestSuite.h" header requires the following parameters:
    - The widths and heights of the two input matrices. To simplify on testing the width of the first matrix is by default the same as the height of the second matrix and is as such inputted only once.
    - The number of tests the program will run.

Test function:
    - NaiveSerialTest: the program will multiply the two matrices using naive algorithm, serial implementation. No padding will be performed.
    - OMPSerialTest: the program will multiply the two matrices using naive algorithm using OMP for parallelization. No padding will be performed.
    - StrassenSerialTest: the program will multiply the two matrices using Strassen algorithm, serial implementation. If either of the two input matrices isn't square and/or max(rowA, colA_rowB, colB) is not a power of 2, the two input matrices will be padded to n x n, with n the smallest power of 2 greater than max(rowA, colA_rowB, colB).
    - StrassenOMP_01Test: the program will multiply the two matrices using Strassen algorithm, parallelizing only the seven sub problems. Padding is the same as StrassenSerialTest.
    - StrassenOMP_02Test: the program will multiply the two matrices using Strassen algorithm, parallelizing only methods for adding, subtractings, splitting and joinning matrices. Padding is the same as StrassenSerialTest.
    - StrassenOMP_03Test the program will multiply the two matrices using Strassen algorithm, parallelizing both the seven subproblems and methods for adding, subtractings, splitting and joinning matrices. Padding is the same as StrassenSerialTest.

## Todo list:
    - Generate log files for recording test results
    - Interface for setting threads, choosing test functions and inputting parameters.