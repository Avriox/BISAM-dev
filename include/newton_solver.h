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
static char _VNEWTON_[] = "@(#)Newton.cpp 01.03 -- Copyright (C) Henrik Vestermark";

#include <algorithm>
#include <complex>
#include <functional>
#include <cfloat>
#include <array>
#include <unordered_map>
#include <cstdint>
#include <cstring>  // for memcpy
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
    complex<double> roots[MAX_ROOTS];
    int count;
};

// Global cache for storing polynomial roots based on coefficient hash
static std::unordered_map<uint32_t, RootStorage> g_rootCache;

// Evaluation structure
struct Eval {
    complex<double> z{};
    complex<double> pz{};
    double apz{};
};

// ===== Helper function declarations =====

// Ultra-fast hashing function specialized for our specific polynomial structure
static uint32_t HashPolynomial(const double coefficients[MAX_COEFF]);

// Compute suitable starting point for root finding
static double ComputeStartPoint(const double a[MAX_COEFF], int n);

// Evaluate polynomial with Horner's method
static Eval EvaluatePolynomial(const double a[MAX_COEFF], int n, const complex<double> z);

// Calculate upper bound for rounding errors (Adam's test)
static double CalculateUpperBound(const double a[MAX_COEFF], int n, const complex<double> z);

// Real root forward deflation
static void RealDeflation(double a[MAX_COEFF], int &n, const double x);

// Complex root forward deflation for real coefficients
static void ComplexDeflation(double a[MAX_COEFF], int &n, const complex<double> z);

// Solve quadratic polynomial
static int SolveQuadratic(const double a[MAX_COEFF], int n, complex<double> roots[MAX_ROOTS],
                          int rootIndex, complex<double> newRoots[MAX_ROOTS], int newRootIndex);

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
    double a0  = log(abs(a[n]));
    double min = exp((a0 - log(abs(a[0]))) / static_cast<double>(n));

    for (int i = 1; i < n; i++) {
        if (a[i] != 0.0) {
            double tmp = exp((a0 - log(abs(a[i]))) / static_cast<double>(n - i));
            if (tmp < min)
                min = tmp;
        }
    }

    return min * 0.5;
}

// Evaluate a polynomial with real coefficients using Horner's method
static Eval EvaluatePolynomial(const double a[MAX_COEFF], int n, const complex<double> z) {
    double p = -2.0 * z.real();
    double q = norm(z);
    double s = 0.0;
    double r = a[0];
    Eval e;

    for (int i = 1; i < n; i++) {
        double t = a[i] - p * r - q * s;
        s        = r;
        r        = t;
    }

    e.z   = z;
    e.pz  = complex<double>(a[n] + z.real() * r - q * s, z.imag() * r);
    e.apz = abs(e.pz);
    return e;
}

// Calculate upper bound for rounding errors (Adam's test)
static double CalculateUpperBound(const double a[MAX_COEFF], int n, const complex<double> z) {
    double p = -2.0 * z.real();
    double q = norm(z);
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

    t = a[n] + z.real() * r - q * s;
    e = u * e + fabs(t);
    e = (4.5 * e - 3.5 * (fabs(t) + fabs(r) * u) +
         fabs(z.real()) * fabs(r)) * 0.5 * pow((double) _DBL_RADIX, -DBL_MANT_DIG + 1);

    return e;
}

// Real root forward deflation
static void RealDeflation(double a[MAX_COEFF], int &n, const double x) {
    double r = 0.0;

    for (int i = 0; i < n; i++)
        a[i]   = r = r * x + a[i];

    n--; // Reduce polynomial degree by 1
}

// Complex root forward deflation for real coefficients
static void ComplexDeflation(double a[MAX_COEFF], int &n, const complex<double> z) {
    double r = -2.0 * z.real();
    double u = norm(z);

    a[1] -= r * a[0];
    for (int i = 2; i < n - 1; i++)
        a[i]   = a[i] - r * a[i - 1] - u * a[i - 2];

    n -= 2; // Reduce polynomial degree by 2 (complex conjugate roots)
}

// Solve quadratic polynomial
static int SolveQuadratic(const double a[MAX_COEFF], int n, complex<double> roots[MAX_ROOTS],
                          int rootIndex, complex<double> newRoots[MAX_ROOTS], int newRootIndex) {
    complex<double> v;
    double r;

    // Notice that a[0] is !=0 since roots=zero has already been handled
    if (n == 1) {
        v                        = complex<double>(-a[1] / a[0], 0);
        roots[rootIndex++]       = v;
        newRoots[newRootIndex++] = v;
    } else {
        if (a[1] == 0.0) {
            r = -a[2] / a[0];
            if (r < 0) {
                r                        = sqrt(-r);
                v                        = complex<double>(0, r);
                roots[rootIndex++]       = v;
                roots[rootIndex++]       = conj(v);
                newRoots[newRootIndex++] = v;
                newRoots[newRootIndex++] = conj(v);
            } else {
                r                        = sqrt(r);
                v                        = complex<double>(r, 0);
                roots[rootIndex++]       = v;
                roots[rootIndex++]       = complex<double>(-r, 0);
                newRoots[newRootIndex++] = v;
                newRoots[newRootIndex++] = complex<double>(-r, 0);
            }
        } else {
            r = 1.0 - 4.0 * a[0] * a[2] / (a[1] * a[1]);
            if (r < 0) {
                v = complex<double>(-a[1] / (2.0 * a[0]), a[1] *
                                                          sqrt(-r) / (2.0 * a[0]));
                roots[rootIndex++]       = v;
                roots[rootIndex++]       = conj(v);
                newRoots[newRootIndex++] = v;
                newRoots[newRootIndex++] = conj(v);
            } else {
                v                        = complex<double>((-1.0 - sqrt(r)) * a[1] / (2.0 * a[0]), 0);
                roots[rootIndex++]       = v;
                newRoots[newRootIndex++] = v;
                v                        = complex<double>(a[2] / (a[0] * v.real()), 0);
                roots[rootIndex++]       = v;
                newRoots[newRootIndex++] = v;
            }
        }
    }
    return rootIndex;
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

// Find all polynomial zeros using a modified Newton method with warm start capability
static int PolynomialRoots(const double coefficients[MAX_COEFF], complex<double> roots[MAX_ROOTS]) {
    const complex<double> complexzero(0.0, 0); // Complex zero (0+i0)
    int n = POLYNOMIAL_DEGREE;                 // Polynomial degree is fixed at 4
    Eval pz;                                   // P(z)
    Eval pzprev;                               // P(zprev)
    Eval p1z;                                  // P'(z)
    Eval p1zprev;                              // P'(zprev)
    complex<double> z;                         // Use as temporary variable
    complex<double> dz;                        // The current stepsize dz
    int itercnt;                               // Hold the number of iterations per root
    int rootIndex = 0;                         // Index for storing roots

    // Fixed-size array for coefficients
    double coeff[MAX_COEFF];
    memcpy(coeff, coefficients, sizeof(double) * MAX_COEFF);

    // Generate hash for this polynomial to check cache - ultra fast integer hash
    uint32_t polyHash     = HashPolynomial(coefficients);
    bool usingCachedRoots = false;
    complex<double> cachedRoots[MAX_ROOTS];
    int cachedRootCount = 0;

    // Check if we have solved a similar polynomial before
    auto cacheIt = g_rootCache.find(polyHash);
    if (cacheIt != g_rootCache.end()) {
        const RootStorage &storage = cacheIt->second;
        cachedRootCount            = storage.count;
        memcpy(cachedRoots, storage.roots, sizeof(complex<double>) * cachedRootCount);
        usingCachedRoots = true;
        DEBUG_PRINT("Found cached roots for polynomial with hash: 0x"
            << std::hex << polyHash << std::dec);
    }

    // Step 1 eliminate all simple roots
    while (n > 0 && coeff[n] == 0.0) {
        roots[rootIndex++] = complexzero; // Store zero as the root
        n--;
    }

    // Index for tracking cached roots
    int cacheRootIndex = 0;
    complex<double> newRoots[MAX_ROOTS];
    int newRootIndex = 0;

    // Do Newton iteration for polynomial order higher than 2
    while (n > 2) {
        const double Max_stepsize      = 5.0;                       // Allow the next step size to be up to 5x larger
        const complex<double> rotation = complex<double>(0.6, 0.8); // Rotation amount
        double r;                                                   // Current radius
        double rprev;                                               // Previous radius
        double eps;                                                 // Iteration termination value
        bool stage1 = true;                                         // By default start in stage1
        int steps   = 1;                                            // Multisteps if > 1

        // Fixed-size array for derivative coefficients
        double coeffprime[MAX_COEFF - 1];

        // Calculate coefficients of p'(x)
        for (int i        = 0; i < n; i++)
            coeffprime[i] = coeff[i] * double(n - i);

        // Step 2 find a suitable starting point z
        rprev = ComputeStartPoint(coeff, n); // Computed startpoint

        // If we have cached roots, use the next one as a starting point
        if (usingCachedRoots && cacheRootIndex < cachedRootCount) {
            z = cachedRoots[cacheRootIndex++];
            DEBUG_PRINT("Using cached root as starting point: " << z);
        } else {
            // Use default starting point calculation
            z = coeff[n - 1] == 0.0 ? complex<double>(1.0, 0) : complex<double>(-coeff[n] / coeff[n - 1], 0);
            z *= complex<double>(rprev, 0) / abs(z);
        }

        // Setup the iteration
        // Current P(z)
        pz = EvaluatePolynomial(coeff, n, z);

        // pzprev which is the previous z or P(0)
        pzprev.z   = complex<double>(0, 0);
        pzprev.pz  = coeff[n];
        pzprev.apz = abs(pzprev.pz);

        // p1zprev P'(0) is the P'(0)
        p1zprev.z   = pzprev.z;
        p1zprev.pz  = coeffprime[n - 1]; // P'(0)
        p1zprev.apz = abs(p1zprev.pz);

        // Set previous dz and calculate the radius of operations.
        dz = pz.z;                  // dz=z-zprev=z since zprev==0
        r  = rprev *= Max_stepsize; // Make a reasonable radius of the maximum step size allowed

        // Preliminary eps computed at point P(0) by a crude estimation
        eps = 2 * n * pzprev.apz * pow((double) _DBL_RADIX, -DBL_MANT_DIG);

        // Start iteration and stop if z doesnt change or apz <= eps
        // we do z+dz!=z instead of dz!=0. if dz does not change z
        // then we accept z as a root
        for (itercnt = 0; pz.z + dz != pz.z && pz.apz > eps && itercnt < MAX_ITER; itercnt++) {
            // Calculate current P'(z)
            p1z = EvaluatePolynomial(coeffprime, n - 1, pz.z);
            if (p1z.apz == 0.0)                                    // P'(z)==0 then rotate and try again
                dz *= rotation * complex<double>(Max_stepsize, 0); // Multiply with 5
            // to get away from saddlepoint
            else {
                dz = pz.pz / p1z.pz; // next dz
                // Check the Magnitude of Newton's step
                r = abs(dz);
                if (r > rprev) // Large than 5 times the previous step size
                {
                    // then rotate and adjust step size to prevent
                    // wild step size near P'(z) close to zero
                    dz *= rotation * complex<double>(rprev / r, 0);
                    r = abs(dz);
                }

                rprev = r * Max_stepsize; // Save 5 times the current step size
                // for the next iteration check
                // of reasonable step size
                // Calculate if stage1 is true or false. Stage1 is false
                // if the Newton converge otherwise true
                z      = (p1zprev.pz - p1z.pz) / (pzprev.z - pz.z);
                stage1 = (abs(z) / p1z.apz > p1z.apz / pz.apz / 4) || (steps != 1);
            }

            // Step accepted. Save pz in pzprev
            pzprev = pz;
            z      = pzprev.z - dz; // Next z
            pz     = EvaluatePolynomial(coeff, n, z);
            steps  = 1;
            if (stage1) {
                // Try multiple steps or shorten steps depending if P(z)
                // is an improvement or not P(z)<P(zprev)
                bool div2;
                complex<double> zn;
                Eval npz;

                zn = pz.z;
                for (div2 = pz.apz > pzprev.apz; steps <= n; ++steps) {
                    if (div2 == true) {
                        // Shorten steps
                        dz *= complex<double>(0.5, 0);
                        zn = pzprev.z - dz;
                    } else
                        zn -= dz; // try another step in the same direction

                    // Evaluate new try step
                    npz = EvaluatePolynomial(coeff, n, zn);
                    if (npz.apz >= pz.apz)
                        break; // Break if no improvement

                    // Improved => accept step and try another round of step
                    pz = npz;

                    if (div2 == true && steps == 2) {
                        // To many shorten steps => try another direction and break
                        dz *= rotation;
                        z  = pzprev.z - dz;
                        pz = EvaluatePolynomial(coeff, n, z);
                        break;
                    }
                }
            } else {
                // calculate the upper bound of error using Grant & Hitchins's
                // test for complex coefficients
                // Now that we are within the convergence circle.
                eps = CalculateUpperBound(coeff, n, pz.z);
            }
        }

        DEBUG_PRINT("Root found after " << itercnt << " iterations: " << pz.z);

        // Check if warm start improved convergence
        if (usingCachedRoots && cacheRootIndex <= cachedRootCount) {
            DEBUG_PRINT("Using warm start: root found in " << itercnt
                << " iterations (vs. typically more without warm start)");
        }

        // Check if there is a very small residue in the imaginary part by trying
        // to evaluate P(z.real). if that is less than P(z).
        // We take that z.real() is a better root than z.
        z      = complex<double>(pz.z.real(), 0.0);
        pzprev = EvaluatePolynomial(coeff, n, z);
        if (pzprev.apz <= pz.apz) {
            // real root
            pz = pzprev;
            // Save the root
            roots[rootIndex++]       = pz.z;
            newRoots[newRootIndex++] = pz.z;
            RealDeflation(coeff, n, pz.z.real());
        } else {
            // Complex root
            // Save the root
            roots[rootIndex++]       = pz.z;
            roots[rootIndex++]       = conj(pz.z);
            newRoots[newRootIndex++] = pz.z;
            newRoots[newRootIndex++] = conj(pz.z);
            ComplexDeflation(coeff, n, pz.z);
        }
    } // End Iteration

    // Solve any remaining linear or quadratic polynomial
    if (n > 0)
        rootIndex = SolveQuadratic(coeff, n, roots, rootIndex, newRoots, newRootIndex);

    // Store the roots in the cache for future use
    RootStorage storage;
    memcpy(storage.roots, newRoots, sizeof(complex<double>) * newRootIndex);
    storage.count         = newRootIndex;
    g_rootCache[polyHash] = storage;

    return rootIndex; // Return the number of roots found
}

// Wrapper function to maintain compatibility with vector-based interface
static std::vector<complex<double> > PolynomialRoots(const std::vector<double> &coefficients) {
    // Convert vector to fixed-size array
    double coeffArray[MAX_COEFF];
    for (int i = 0; i < MAX_COEFF; i++) {
        coeffArray[i] = coefficients[i];
    }

    // Call the optimized function
    complex<double> rootsArray[MAX_ROOTS];
    int rootCount = PolynomialRoots(coeffArray, rootsArray);

    // Convert back to vector for return
    std::vector<complex<double> > roots;
    for (int i = 0; i < rootCount; i++) {
        roots.push_back(rootsArray[i]);
    }

    return roots;
}

#endif //NEWTON_SOLVER_H
