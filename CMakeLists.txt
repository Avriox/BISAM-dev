cmake_minimum_required(VERSION 3.30)
project(BISAM)

set(CMAKE_CXX_STANDARD 17)

# Defines DEV - useful for conditional compilation during development
add_compile_definitions(DEV)

# Special settings to show full command line during build (optional)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native -flto -funroll-loops -fomit-frame-pointer")

# macOS-specific compiler approach
if (APPLE)
    # First try to find OpenMP normally
    find_package(OpenMP QUIET)

    # If not found, use the macOS-specific approach
    if (NOT OpenMP_FOUND)
        # Check if Homebrew's libomp is installed
        if (EXISTS "/opt/homebrew/opt/libomp")
            # For macOS, we add the flags directly rather than using OpenMP_* variables
            add_compile_options(-Xpreprocessor -fopenmp)
            include_directories("/opt/homebrew/opt/libomp/include")
            link_directories("/opt/homebrew/opt/libomp/lib")
            set(OPENMP_ENABLED TRUE)
            message(STATUS "OpenMP enabled with Homebrew libomp")
        else ()
            message(STATUS "OpenMP not found - building without parallel support")
        endif ()
    else ()
        message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
        set(OPENMP_ENABLED TRUE)
    endif ()
else ()
    # For non-macOS systems, use standard OpenMP detection
    find_package(OpenMP)
    if (OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
        set(OPENMP_ENABLED TRUE)
    else ()
        message(STATUS "OpenMP not found - building without parallel support")
    endif ()
endif ()

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

# Find OpenBLAS (or other BLAS/LAPACK implementations)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Print config details for verification
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

# Define the executable and list all source files
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
        include/globals.h
        include/fast_polynomial_solver.h
        src/fast_polynomial_solver.cpp
)

# Include all necessary directories for header lookup
target_include_directories(BISAM PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/lib/lasso
        ${R_INCLUDE_DIR}
        ${RCPPARMADILLO_INCLUDE_DIR}
        ${RCPP_INCLUDE_DIR}
        ${MOMBF_SRC_DIR}
)

# Link against required libraries
target_link_libraries(BISAM PRIVATE
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        -L${R_LIB_DIR} -lR
)

# Add OpenMP-related settings
if (OPENMP_ENABLED)
    # If using standard OpenMP detection
    if (OpenMP_CXX_FOUND)
        target_compile_options(BISAM PRIVATE ${OpenMP_CXX_FLAGS})
        target_link_libraries(BISAM PRIVATE ${OpenMP_CXX_LIBRARIES})
        # If using macOS direct approach
    elseif (APPLE)
        target_link_libraries(BISAM PRIVATE -lomp)
    endif ()

    # Add a definition so the C++ code can check if OpenMP is enabled
    target_compile_definitions(BISAM PRIVATE HAVE_OPENMP)
endif ()