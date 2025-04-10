cmake_minimum_required(VERSION 3.30)
project(BISAM)

set(CMAKE_CXX_STANDARD 17)

add_compile_definitions(DEV)

# Special settings to show full command line
set(CMAKE_VERBOSE_MAKEFILE ON)

# Find Rscript (more suitable for scripting than R)
find_program(RSCRIPT_EXECUTABLE Rscript)
if (NOT RSCRIPT_EXECUTABLE)
    message(FATAL_ERROR "Rscript executable not found")
endif ()

# Get R include directory
execute_process(
        COMMAND ${RSCRIPT_EXECUTABLE} -e "cat(R.home('include'))"
        OUTPUT_VARIABLE R_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get R library directory
execute_process(
        COMMAND ${RSCRIPT_EXECUTABLE} -e "cat(R.home('lib'))"
        OUTPUT_VARIABLE R_LIB_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get RcppArmadillo include directory
execute_process(
        COMMAND ${RSCRIPT_EXECUTABLE} -e "cat(system.file('include', package='RcppArmadillo'))"
        OUTPUT_VARIABLE RCPPARMADILLO_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get Rcpp include directory
execute_process(
        COMMAND ${RSCRIPT_EXECUTABLE} -e "cat(system.file('include', package='Rcpp'))"
        OUTPUT_VARIABLE RCPP_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Find OpenBLAS (try this before finding Armadillo)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Find regular Armadillo (might still be needed)
find_package(Armadillo REQUIRED)

# Print config details
message(STATUS "R include: ${R_INCLUDE_DIR}")
message(STATUS "R lib: ${R_LIB_DIR}")
message(STATUS "RcppArmadillo include: ${RCPPARMADILLO_INCLUDE_DIR}")
message(STATUS "Rcpp include: ${RCPP_INCLUDE_DIR}")
message(STATUS "Armadillo includes: ${ARMADILLO_INCLUDE_DIRS}")
message(STATUS "Armadillo libs: ${ARMADILLO_LIBRARIES}")
message(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")

# Define the directory containing the mombf files
set(MOMBF_SRC_DIR mombf/src)

# Find all .h and .cpp files in the mombf source directory
file(GLOB_RECURSE MOMBF_SOURCES
        "${MOMBF_SRC_DIR}/*.h"
        "${MOMBF_SRC_DIR}/*.cpp"
)

add_executable(BISAM main.cpp
        src/controlled-simulation.h
        src/mombf-bridge.h
        ${MOMBF_SOURCES}
        lib/lasso/LassoRegression.h
        lib/lasso/LassoRegression.cpp
        lib/lasso/matrix.cpp
        lib/lasso/matrix.h
        src/model-selection-strategy.h
        src/exec_timer.h
        src/b_ism.h
        src/bisam_rcpp.cpp
        src/b_ism.cpp
        src/bisam_rcpp.h
)

# Include all necessary directories
target_include_directories(BISAM PRIVATE
        ${R_INCLUDE_DIR}
        ${RCPPARMADILLO_INCLUDE_DIR}
        ${RCPP_INCLUDE_DIR}
        ${ARMADILLO_INCLUDE_DIRS}
        ${MOMBF_SRC_DIR} # Include the mombf source directory for header files
)

# Link against libraries
target_link_libraries(BISAM PRIVATE
        ${ARMADILLO_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        -L${R_LIB_DIR} -lR
)