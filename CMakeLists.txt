cmake_minimum_required(VERSION 3.30)
project(BISAM)

set(CMAKE_CXX_STANDARD 17)

# Defines DEV
add_compile_definitions(DEV)

# Find OpenMP
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
else ()
    message(STATUS "OpenMP not found - building without parallel support")
endif ()

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


# Print config details
message(STATUS "R include: ${R_INCLUDE_DIR}")
message(STATUS "R lib: ${R_LIB_DIR}")
message(STATUS "RcppArmadillo include: ${RCPPARMADILLO_INCLUDE_DIR}")
message(STATUS "Rcpp include: ${RCPP_INCLUDE_DIR}")
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
        ${MOMBF_SOURCES}
        lib/lasso/LassoRegression.h
        lib/lasso/LassoRegression.cpp
        lib/lasso/matrix.cpp
        lib/lasso/matrix.h
        src/bisam_rcpp.cpp
        include/bisam_types.h
        include/biasm_model.h
        include/mombf_bridge.h
        include/modelselection_strategy.h
        include/utils.h
        src/biasm_model.cpp
        src/mombf_bridge.cpp
        src/modelselection_strategy.cpp
        src/utils.cpp
)

# Include all necessary directories
target_include_directories(BISAM PRIVATE
        ${R_INCLUDE_DIR}
        ${RCPPARMADILLO_INCLUDE_DIR}
        ${RCPP_INCLUDE_DIR}
        ${MOMBF_SRC_DIR} # Include the mombf source directory for header files
)

# Link against libraries
target_link_libraries(BISAM PRIVATE
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        -L${R_LIB_DIR} -lR
)

# Add OpenMP flags if found
if (OpenMP_CXX_FOUND)
    target_compile_options(BISAM PRIVATE ${OpenMP_CXX_FLAGS})
    target_link_libraries(BISAM PRIVATE ${OpenMP_CXX_LIBRARIES})
    # Add a definition so the code knows OpenMP is available
    target_compile_definitions(BISAM PRIVATE HAVE_OPENMP)
endif ()