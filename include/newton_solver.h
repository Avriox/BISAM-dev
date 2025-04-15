#ifndef NEWTON_SOLVER_H
#define NEWTON_SOLVER_H

/*
 *******************************************************************************
 *
 *                       Copyright (c) 2023
 *                       Henrik Vestermark
 *                       Denmark, USA
 *
 *                       All Rights Reserved
 *
 *   This source file is subject to the terms and conditions of
 *   Henrik Vestermark Software License Agreement which restricts the manner
 *   in which it may be used.
 *
 *******************************************************************************
*/
/*
 *******************************************************************************
 *
 * Module name     :   Newton.cpp
 * Module ID Nbr   :
 * Description     :   Solve n degree polynomial using Newton's (Madsen) method
 * --------------------------------------------------------------------------
 * Change Record   :
 *
 * Version      Author/Date     Description of changes
 * -------     -------------  ----------------------
 * 01.01        HVE/24Sep2023  Initial release
 * 01.02        ---/15Apr2025  Added warm start capability
 * 01.03        ---/15Apr2025  Optimized with fixed-size arrays
 * 01.04        ---/15Apr2025  Direct C-style interface with no vector operations
 * 01.05        ---/15Apr2025  Eliminated std::complex operations
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/

// Define to enable debug printing. Comment out to disable.
// #define DEBUG_PRINTING

#ifdef DEBUG_PRINTING
#include <iostream>
#define DEBUG_PRINT(x) std::cout << x << std::endl
#else
#define DEBUG_PRINT(x)
#endif

// define version string
static char _VNEWTON_[] = "@(#)Newton.cpp 01.05 -- Copyright (C) Henrik Vestermark";

#include <algorithm>
#include <functional>
#include <cfloat>
#include <unordered_map>
#include <cstdint>
#include <cstring>  // for memcpy
#include <cmath>    // for sqrt, abs, etc.
#include "PolynomialRootFinder.h"

//[BISAM] _DBL_RADIX was not defined
# define _DBL_RADIX FLT_RADIX

using namespace std;
constexpr int MAX_ITER          = 50;
constexpr int POLYNOMIAL_DEGREE = 4;                     // We always have degree 4 polynomials
constexpr int MAX_COEFF         = POLYNOMIAL_DEGREE + 1; // Number of coefficients in a degree 4 polynomial
constexpr int MAX_ROOTS         = POLYNOMIAL_DEGREE;     // Maximum number of roots

// ===== Global cache and helper structs =====

// Structure to store roots with a count of valid roots
struct RootStorage {
    double real[MAX_ROOTS]; // Real parts of the roots
    double imag[MAX_ROOTS]; // Imaginary parts of the roots
    int count;              // Number of valid roots
};

// Global cache for storing polynomial roots based on coefficient hash
static std::unordered_map<uint32_t, RootStorage> g_rootCache;

// Structure representing a complex number without using std::complex
struct Complex {
    double re; // Real part
    double im; // Imaginary part

    // Fast norm (squared magnitude)
    inline double norm() const {
        return re * re + im * im;
    }

    // Fast absolute value (magnitude)
    inline double abs() const {
        return sqrt(norm());
    }
};

// Structure for polynomial evaluation
struct Eval {
    double z_re;  // Real part of z
    double z_im;  // Imaginary part of z
    double pz_re; // Real part of P(z)
    double pz_im; // Imaginary part of P(z)
    double apz;   // |P(z)| - absolute value
};

// ===== Complex arithmetic operations =====

// Complex number addition
inline void complex_add(double &result_re, double &result_im,
                        double a_re, double a_im,
                        double b_re, double b_im) {
    result_re = a_re + b_re;
    result_im = a_im + b_im;
}

// Complex number subtraction
inline void complex_sub(double &result_re, double &result_im,
                        double a_re, double a_im,
                        double b_re, double b_im) {
    result_re = a_re - b_re;
    result_im = a_im - b_im;
}

// Complex number multiplication
inline void complex_mul(double &result_re, double &result_im,
                        double a_re, double a_im,
                        double b_re, double b_im) {
    double tmp_re = a_re * b_re - a_im * b_im;
    double tmp_im = a_re * b_im + a_im * b_re;
    result_re     = tmp_re;
    result_im     = tmp_im;
}

// Complex number division
inline void complex_div(double &result_re, double &result_im,
                        double a_re, double a_im,
                        double b_re, double b_im) {
    double denom  = b_re * b_re + b_im * b_im;
    double tmp_re = (a_re * b_re + a_im * b_im) / denom;
    double tmp_im = (a_im * b_re - a_re * b_im) / denom;
    result_re     = tmp_re;
    result_im     = tmp_im;
}

// Complex number absolute value (magnitude)
inline double complex_abs(double re, double im) {
    return sqrt(re * re + im * im);
}

// Complex number conjugate
inline void complex_conj(double &result_re, double &result_im,
                         double a_re, double a_im) {
    result_re = a_re;
    result_im = -a_im;
}

// Scale complex number by a real value
inline void complex_scale(double &result_re, double &result_im,
                          double a_re, double a_im,
                          double scale) {
    result_re = a_re * scale;
    result_im = a_im * scale;
}

// ===== Helper function declarations =====

// Ultra-fast hashing function specialized for our specific polynomial structure
static uint32_t HashPolynomial(const double coefficients[MAX_COEFF]);

// Compute suitable starting point for root finding
static double ComputeStartPoint(const double a[MAX_COEFF], int n);

// Evaluate polynomial with Horner's method
static Eval EvaluatePolynomial(const double a[MAX_COEFF], int n, double z_re, double z_im);

// Calculate upper bound for rounding errors (Adam's test)
static double CalculateUpperBound(const double a[MAX_COEFF], int n, double z_re, double z_im);

// Real root forward deflation
static void RealDeflation(double a[MAX_COEFF], int &n, double x);

// Complex root forward deflation for real coefficients
static void ComplexDeflation(double a[MAX_COEFF], int &n, double z_re, double z_im);

// Solve quadratic polynomial
static int SolveQuadratic(const double a[MAX_COEFF], int n,
                          double real_roots[MAX_ROOTS], double imag_roots[MAX_ROOTS],
                          int rootIndex, double new_real[MAX_ROOTS], double new_imag[MAX_ROOTS],
                          int newRootIndex);

// Root cache management functions
static void ClearRootCache();

static void PrintRootCacheStats();

static void ReserveRootCache(size_t expectedSize);

// ===== Implementation of helper functions =====

// Ultra-fast hashing function specialized for our specific polynomial structure
static uint32_t HashPolynomial(const double coefficients[MAX_COEFF]) {
    // Get the integer parts of coefficients at positions 4, 1, and 0
    int c4_int = static_cast<int>(coefficients[4]);
    int c1_int = static_cast<int>(coefficients[1]);
    int c0_int = static_cast<int>(coefficients[0]);

    // Convert to positive values for consistent hashing
    uint32_t c4 = static_cast<uint32_t>(c4_int < 0 ? -c4_int : c4_int);
    uint32_t c1 = static_cast<uint32_t>(c1_int < 0 ? -c1_int : c1_int);
    uint32_t c0 = static_cast<uint32_t>(c0_int < 0 ? -c0_int : c0_int);

    // Create hash by bit-packing (very fast)
    // We reserve 10 bits for each value (allows values up to 1023)
    // Add sign bits in highest positions
    uint32_t hash = (c0 << 20) | (c1 << 10) | c4;
    // Include sign information in the two highest bits
    if (c0_int < 0) hash |= (1U << 30);
    if (c1_int < 0) hash |= (1U << 31);

    return hash;
}

// Compute the next starting point based on the polynomial coefficients
static double ComputeStartPoint(const double a[MAX_COEFF], int n) {
    double a0  = log(fabs(a[n]));
    double min = exp((a0 - log(fabs(a[0]))) / static_cast<double>(n));

    for (int i = 1; i < n; i++) {
        if (a[i] != 0.0) {
            double tmp = exp((a0 - log(fabs(a[i]))) / static_cast<double>(n - i));
            if (tmp < min)
                min = tmp;
        }
    }

    return min * 0.5;
}

// Evaluate a polynomial with real coefficients using Horner's method
// Directly calculates real and imaginary parts without using complex numbers
static Eval EvaluatePolynomial(const double a[MAX_COEFF], int n, double z_re, double z_im) {
    double p = -2.0 * z_re;
    double q = z_re * z_re + z_im * z_im; // Norm of z
    double s = 0.0;
    double r = a[0];
    Eval e;

    for (int i = 1; i < n; i++) {
        double t = a[i] - p * r - q * s;
        s        = r;
        r        = t;
    }

    // Calculate P(z)
    e.z_re  = z_re;
    e.z_im  = z_im;
    e.pz_re = a[n] + z_re * r - q * s;
    e.pz_im = z_im * r;
    e.apz   = sqrt(e.pz_re * e.pz_re + e.pz_im * e.pz_im); // |P(z)|

    return e;
}

// Calculate upper bound for rounding errors (Adam's test)
static double CalculateUpperBound(const double a[MAX_COEFF], int n, double z_re, double z_im) {
    double p = -2.0 * z_re;
    double q = z_re * z_re + z_im * z_im; // Norm of z
    double u = sqrt(q);
    double s = 0.0;
    double r = a[0];
    double e = fabs(r) * (3.5 / 4.5);
    double t;

    for (int i = 1; i < n; i++) {
        t = a[i] - p * r - q * s;
        s = r;
        r = t;
        e = u * e + fabs(t);
    }

    t = a[n] + z_re * r - q * s;
    e = u * e + fabs(t);
    e = (4.5 * e - 3.5 * (fabs(t) + fabs(r) * u) +
         fabs(z_re) * fabs(r)) * 0.5 * pow((double) _DBL_RADIX, -DBL_MANT_DIG + 1);

    return e;
}

// Real root forward deflation
static void RealDeflation(double a[MAX_COEFF], int &n, double x) {
    double r = 0.0;

    for (int i = 0; i < n; i++)
        a[i]   = r = r * x + a[i];

    n--; // Reduce polynomial degree by 1
}

// Complex root forward deflation for real coefficients
static void ComplexDeflation(double a[MAX_COEFF], int &n, double z_re, double z_im) {
    double r = -2.0 * z_re;
    double u = z_re * z_re + z_im * z_im; // Norm of z

    a[1] -= r * a[0];
    for (int i = 2; i < n - 1; i++)
        a[i]   = a[i] - r * a[i - 1] - u * a[i - 2];

    n -= 2; // Reduce polynomial degree by 2 (complex conjugate roots)
}

// Solve quadratic polynomial - directly write real and imaginary parts
static int SolveQuadratic(const double a[MAX_COEFF], int n,
                          double real_roots[MAX_ROOTS], double imag_roots[MAX_ROOTS],
                          int rootIndex, double new_real[MAX_ROOTS], double new_imag[MAX_ROOTS],
                          int newRootIndex) {
    double real_part, imag_part;
    double r;

    // Notice that a[0] is !=0 since roots=zero has already been handled
    if (n == 1) {
        real_part = -a[1] / a[0];
        imag_part = 0.0;

        real_roots[rootIndex]  = real_part;
        imag_roots[rootIndex]  = imag_part;
        new_real[newRootIndex] = real_part;
        new_imag[newRootIndex] = imag_part;

        return rootIndex + 1;
    } else {
        if (a[1] == 0.0) {
            r = -a[2] / a[0];
            if (r < 0) {
                r         = sqrt(-r);
                real_part = 0.0;
                imag_part = r;

                real_roots[rootIndex]     = real_part;
                imag_roots[rootIndex]     = imag_part;
                real_roots[rootIndex + 1] = real_part;
                imag_roots[rootIndex + 1] = -imag_part;

                new_real[newRootIndex]     = real_part;
                new_imag[newRootIndex]     = imag_part;
                new_real[newRootIndex + 1] = real_part;
                new_imag[newRootIndex + 1] = -imag_part;

                return rootIndex + 2;
            } else {
                r = sqrt(r);

                real_roots[rootIndex]     = r;
                imag_roots[rootIndex]     = 0.0;
                real_roots[rootIndex + 1] = -r;
                imag_roots[rootIndex + 1] = 0.0;

                new_real[newRootIndex]     = r;
                new_imag[newRootIndex]     = 0.0;
                new_real[newRootIndex + 1] = -r;
                new_imag[newRootIndex + 1] = 0.0;

                return rootIndex + 2;
            }
        } else {
            r = 1.0 - 4.0 * a[0] * a[2] / (a[1] * a[1]);
            if (r < 0) {
                real_part = -a[1] / (2.0 * a[0]);
                imag_part = a[1] * sqrt(-r) / (2.0 * a[0]);

                real_roots[rootIndex]     = real_part;
                imag_roots[rootIndex]     = imag_part;
                real_roots[rootIndex + 1] = real_part;
                imag_roots[rootIndex + 1] = -imag_part;

                new_real[newRootIndex]     = real_part;
                new_imag[newRootIndex]     = imag_part;
                new_real[newRootIndex + 1] = real_part;
                new_imag[newRootIndex + 1] = -imag_part;

                return rootIndex + 2;
            } else {
                real_part = (-1.0 - sqrt(r)) * a[1] / (2.0 * a[0]);
                imag_part = 0.0;

                real_roots[rootIndex] = real_part;
                imag_roots[rootIndex] = imag_part;

                new_real[newRootIndex] = real_part;
                new_imag[newRootIndex] = imag_part;

                real_part = a[2] / (a[0] * real_part);

                real_roots[rootIndex + 1] = real_part;
                imag_roots[rootIndex + 1] = 0.0;

                new_real[newRootIndex + 1] = real_part;
                new_imag[newRootIndex + 1] = 0.0;

                return rootIndex + 2;
            }
        }
    }
}

// Root cache management functions
static void ClearRootCache() {
    g_rootCache.clear();
    DEBUG_PRINT("Root cache cleared");
}

static void PrintRootCacheStats() {
    DEBUG_PRINT("Root cache contains " << g_rootCache.size() << " entries");
    DEBUG_PRINT("Root cache bucket count: " << g_rootCache.bucket_count());
    DEBUG_PRINT("Root cache load factor: " << g_rootCache.load_factor());
}

static void ReserveRootCache(size_t expectedSize) {
    g_rootCache.reserve(expectedSize);
    DEBUG_PRINT("Reserved cache for " << expectedSize << " polynomials");
}

// ===== Main polynomial roots function =====

// C-style interface that directly writes to real and imaginary arrays
// Optimized version without std::complex operations
static int PolynomialRoots(
    const double *coefficient_vector_ptr,
    int degree,
    double *real_zero_vector_ptr,
    double *imaginary_zero_vector_ptr) {
    int n = degree;      // Polynomial degree
    Eval pz;             // P(z)
    Eval pzprev;         // P(zprev)
    Eval p1z;            // P'(z)
    Eval p1zprev;        // P'(zprev)
    double z_re, z_im;   // Current z (replaces complex<double> z)
    double dz_re, dz_im; // Current step size (replaces complex<double> dz)
    int itercnt;         // Hold the number of iterations per root
    int rootIndex = 0;   // Index for storing roots

    // Constants for complex rotation instead of creating complex numbers
    const double rotation_re = 0.6;
    const double rotation_im = 0.8;

    // Fixed-size array for coefficients
    double coeff[MAX_COEFF];
    memcpy(coeff, coefficient_vector_ptr, sizeof(double) * (degree + 1));

    // Generate hash for this polynomial to check cache - ultra fast integer hash
    uint32_t polyHash     = HashPolynomial(coeff);
    bool usingCachedRoots = false;
    double cachedReal[MAX_ROOTS];
    double cachedImag[MAX_ROOTS];
    int cachedRootCount = 0;

    // Check if we have solved a similar polynomial before
    auto cacheIt = g_rootCache.find(polyHash);
    if (cacheIt != g_rootCache.end()) {
        const RootStorage &storage = cacheIt->second;
        cachedRootCount            = storage.count;
        memcpy(cachedReal, storage.real, sizeof(double) * cachedRootCount);
        memcpy(cachedImag, storage.imag, sizeof(double) * cachedRootCount);
        usingCachedRoots = true;
        DEBUG_PRINT("Found cached roots for polynomial with hash: 0x"
            << std::hex << polyHash << std::dec);
    }

    // Step 1 eliminate all simple roots
    while (n > 0 && coeff[n] == 0.0) {
        real_zero_vector_ptr[rootIndex]      = 0.0; // Store real part of zero root
        imaginary_zero_vector_ptr[rootIndex] = 0.0; // Store imaginary part of zero root
        rootIndex++;
        n--;
    }

    // Index for tracking cached roots
    int cacheRootIndex = 0;
    double newReal[MAX_ROOTS];
    double newImag[MAX_ROOTS];
    int newRootIndex = 0;

    // Do Newton iteration for polynomial order higher than 2
    while (n > 2) {
        const double Max_stepsize = 5.0; // Allow the next step size to be up to 5x larger
        double r;                        // Current radius
        double rprev;                    // Previous radius
        double eps;                      // Iteration termination value
        bool stage1 = true;              // By default start in stage1
        int steps   = 1;                 // Multisteps if > 1

        // Fixed-size array for derivative coefficients
        double coeffprime[MAX_COEFF - 1];

        // Calculate coefficients of p'(x)
        for (int i        = 0; i < n; i++)
            coeffprime[i] = coeff[i] * double(n - i);

        // Step 2 find a suitable starting point z
        rprev = ComputeStartPoint(coeff, n); // Computed startpoint

        // If we have cached roots, use the next one as a starting point
        if (usingCachedRoots && cacheRootIndex < cachedRootCount) {
            z_re = cachedReal[cacheRootIndex];
            z_im = cachedImag[cacheRootIndex];
            cacheRootIndex++;
            DEBUG_PRINT("Using cached root as starting point: " << z_re << " + " << z_im << "i");
        } else {
            // Use default starting point calculation
            if (coeff[n - 1] == 0.0) {
                z_re = 1.0;
                z_im = 0.0;
            } else {
                z_re = -coeff[n] / coeff[n - 1];
                z_im = 0.0;
            }

            // Scale by rprev / |z|
            double abs_z = sqrt(z_re * z_re + z_im * z_im);
            if (abs_z > 0.0) {
                double scale = rprev / abs_z;
                z_re *= scale;
                z_im *= scale;
            }
        }

        // Setup the iteration
        // Current P(z)
        pz = EvaluatePolynomial(coeff, n, z_re, z_im);

        // pzprev which is the previous z or P(0)
        pzprev.z_re  = 0.0;
        pzprev.z_im  = 0.0;
        pzprev.pz_re = coeff[n];
        pzprev.pz_im = 0.0;
        pzprev.apz   = fabs(coeff[n]);

        // p1zprev P'(0) is the P'(0)
        p1zprev.z_re  = pzprev.z_re;
        p1zprev.z_im  = pzprev.z_im;
        p1zprev.pz_re = coeffprime[n - 1];
        p1zprev.pz_im = 0.0;
        p1zprev.apz   = fabs(coeffprime[n - 1]);

        // Set previous dz and calculate the radius of operations.
        dz_re = pz.z_re;
        dz_im = pz.z_im;
        r     = rprev *= Max_stepsize; // Make a reasonable radius of the maximum step size allowed

        // Preliminary eps computed at point P(0) by a crude estimation
        eps = 2 * n * pzprev.apz * pow((double) _DBL_RADIX, -DBL_MANT_DIG);

        // Start iteration
        for (itercnt = 0;
             // Check if z + dz != z (equivalent to dz != 0 in complex arithmetic)
             (pz.z_re + dz_re != pz.z_re || pz.z_im + dz_im != pz.z_im) &&
             pz.apz > eps && itercnt < MAX_ITER;
             itercnt++) {
            // Calculate current P'(z)
            p1z = EvaluatePolynomial(coeffprime, n - 1, pz.z_re, pz.z_im);

            // If P'(z)==0 then rotate and try again
            if (p1z.apz == 0.0) {
                // Multiply dz by rotation * Max_stepsize to get away from saddlepoint
                double tmp_re, tmp_im;
                complex_mul(tmp_re, tmp_im, dz_re, dz_im, rotation_re, rotation_im);
                complex_scale(dz_re, dz_im, tmp_re, tmp_im, Max_stepsize);
            } else {
                // Calculate next dz = pz.pz / p1z.pz using complex division
                complex_div(dz_re, dz_im, pz.pz_re, pz.pz_im, p1z.pz_re, p1z.pz_im);

                // Check the Magnitude of Newton's step
                r = sqrt(dz_re * dz_re + dz_im * dz_im);

                if (r > rprev) {
                    // Too large step - rotate and adjust step size
                    double tmp_re, tmp_im, scaled_re, scaled_im;
                    complex_mul(tmp_re, tmp_im, dz_re, dz_im, rotation_re, rotation_im);
                    complex_scale(scaled_re, scaled_im, tmp_re, tmp_im, rprev / r);
                    dz_re = scaled_re;
                    dz_im = scaled_im;
                    r     = sqrt(dz_re * dz_re + dz_im * dz_im);
                }

                // Save 5 times the current step size for the next iteration check
                rprev = r * Max_stepsize;

                // Calculate if stage1 is true or false
                // Stage1 is false if the Newton converges, otherwise true
                double z_re_tmp, z_im_tmp;
                complex_sub(z_re_tmp, z_im_tmp, p1zprev.pz_re, p1zprev.pz_im, p1z.pz_re, p1z.pz_im);
                double denom_re, denom_im;
                complex_sub(denom_re, denom_im, pzprev.z_re, pzprev.z_im, pz.z_re, pz.z_im);

                // z = (p1zprev.pz - p1z.pz) / (pzprev.z - pz.z)
                complex_div(z_re_tmp, z_im_tmp, z_re_tmp, z_im_tmp, denom_re, denom_im);

                // Calculate: (abs(z) / p1z.apz > p1z.apz / pz.apz / 4)
                double abs_z_tmp = sqrt(z_re_tmp * z_re_tmp + z_im_tmp * z_im_tmp);
                stage1           = (abs_z_tmp / p1z.apz > p1z.apz / pz.apz / 4) || (steps != 1);
            }

            // Step accepted. Save pz in pzprev
            pzprev = pz;

            // Next z = pzprev.z - dz
            complex_sub(z_re, z_im, pzprev.z_re, pzprev.z_im, dz_re, dz_im);

            // Evaluate polynomial at new z
            pz    = EvaluatePolynomial(coeff, n, z_re, z_im);
            steps = 1;

            if (stage1) {
                // Try multiple steps or shorten steps depending if P(z)
                // is an improvement or not P(z)<P(zprev)
                bool div2;
                double zn_re = pz.z_re;
                double zn_im = pz.z_im;
                Eval npz;

                for (div2 = pz.apz > pzprev.apz; steps <= n; ++steps) {
                    if (div2 == true) {
                        // Shorten steps
                        dz_re *= 0.5;
                        dz_im *= 0.5;
                        complex_sub(zn_re, zn_im, pzprev.z_re, pzprev.z_im, dz_re, dz_im);
                    } else {
                        // Try another step in the same direction
                        complex_sub(zn_re, zn_im, zn_re, zn_im, dz_re, dz_im);
                    }

                    // Evaluate new try step
                    npz = EvaluatePolynomial(coeff, n, zn_re, zn_im);
                    if (npz.apz >= pz.apz)
                        break; // Break if no improvement

                    // Improved => accept step and try another round of step
                    pz = npz;

                    if (div2 == true && steps == 2) {
                        // Too many shorten steps => try another direction and break
                        double tmp_re, tmp_im;
                        complex_mul(tmp_re, tmp_im, dz_re, dz_im, rotation_re, rotation_im);
                        dz_re = tmp_re;
                        dz_im = tmp_im;
                        complex_sub(z_re, z_im, pzprev.z_re, pzprev.z_im, dz_re, dz_im);
                        pz = EvaluatePolynomial(coeff, n, z_re, z_im);
                        break;
                    }
                }
            } else {
                // Calculate the upper bound of error using Grant & Hitchins's test
                eps = CalculateUpperBound(coeff, n, pz.z_re, pz.z_im);
            }
        }

        DEBUG_PRINT("Root found after " << itercnt << " iterations: "
            << pz.z_re << " + " << pz.z_im << "i");

        // Check if warm start improved convergence
        if (usingCachedRoots && cacheRootIndex <= cachedRootCount) {
            DEBUG_PRINT("Using warm start: root found in " << itercnt
                << " iterations (vs. typically more without warm start)");
        }

        // Check if there is a very small residue in the imaginary part
        // Try to evaluate P(z.real), if that is less than P(z)
        // Take that z.real() is a better root than z
        Eval pztest = EvaluatePolynomial(coeff, n, pz.z_re, 0.0);

        if (pztest.apz <= pz.apz) {
            // Real root
            pz = pztest;

            // Save the root directly to output arrays
            real_zero_vector_ptr[rootIndex]      = pz.z_re;
            imaginary_zero_vector_ptr[rootIndex] = 0.0;
            rootIndex++;

            // Also save to the cache arrays
            newReal[newRootIndex] = pz.z_re;
            newImag[newRootIndex] = 0.0;
            newRootIndex++;

            RealDeflation(coeff, n, pz.z_re);
        } else {
            // Complex conjugate pair of roots

            // Save both roots directly to output arrays
            real_zero_vector_ptr[rootIndex]      = pz.z_re;
            imaginary_zero_vector_ptr[rootIndex] = pz.z_im;
            rootIndex++;

            real_zero_vector_ptr[rootIndex]      = pz.z_re;
            imaginary_zero_vector_ptr[rootIndex] = -pz.z_im;
            rootIndex++;

            // Also save to the cache arrays
            newReal[newRootIndex] = pz.z_re;
            newImag[newRootIndex] = pz.z_im;
            newRootIndex++;

            newReal[newRootIndex] = pz.z_re;
            newImag[newRootIndex] = -pz.z_im;
            newRootIndex++;

            ComplexDeflation(coeff, n, pz.z_re, pz.z_im);
        }
    } // End Iteration

    // Solve any remaining linear or quadratic polynomial
    if (n > 0)
        rootIndex = SolveQuadratic(coeff, n, real_zero_vector_ptr, imaginary_zero_vector_ptr,
                                   rootIndex, newReal, newImag, newRootIndex);

    // Store the roots in the cache for future use
    RootStorage storage;
    memcpy(storage.real, newReal, sizeof(double) * newRootIndex);
    memcpy(storage.imag, newImag, sizeof(double) * newRootIndex);
    storage.count         = newRootIndex;
    g_rootCache[polyHash] = storage;

    return rootIndex; // Return the number of roots found
}

#endif //NEWTON_SOLVER_H
