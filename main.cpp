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
    TestSuite::StrassenOMP_02Test(4096,4096,4096,5);
    TestSuite::StrassenOMP_01Test(4096,4096,4096,5);
    TestSuite::StrassenOMP_03Test(4096,4096,4096,5);
    return 0;
}