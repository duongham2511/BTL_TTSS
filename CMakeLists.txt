cmake_minimum_required(VERSION 3.0.0)
project(StrassenOMP VERSION 0.1.0)

include(CTest)
enable_testing()

find_package(OpenMP REQUIRED)
add_executable(StrassenOMP main.cpp)
target_link_libraries(StrassenOMP -fopenmp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
