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
    omp_set_num_threads(8);
    TestSuite::StrassenOMP_01Test(8192,8192,8192,3);
    return 0;
}