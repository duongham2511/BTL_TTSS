#include <omp.h>
#include <stdio.h>
#include <math.h>
#include "Utils.h"
#include "MatOp.h"
#include "MatOp_P.h"
#include "Strassen.h"
#include "TestSuite.h"


int main(int argc, char **argv)
{
    omp_set_num_threads(32);
    TestSuite::StrassenOMP_02Test(2048,2048,2048,3);
    TestSuite::StrassenOMP_01Test(2048,2048,2048,3);
    return 0;
}