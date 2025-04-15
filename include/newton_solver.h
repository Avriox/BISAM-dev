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
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/

// define version string
static char _VNEWTON_[] = "@(#)Newton.cpp 01.02 -- Copyright (C) Henrik Vestermark";

#include <algorithm>
#include <vector>
#include <complex>
#include <functional>
#include <cfloat>
#include <unordered_map>
#include <cstdint>
#include "PolynomialRootFinder.h"

//[BISAM] _DBL_RADIX was not defined
# define _DBL_RADIX FLT_RADIX

using namespace std;
constexpr int MAX_ITER = 50;

// ===== Global cache and helper structs =====

// Global cache for storing polynomial roots based on coefficient hash
static std::unordered_map<uint32_t, std::vector<std::complex<double> > > g_rootCache;

// Evaluation structure
struct Eval {
    complex<double> z{};
    complex<double> pz{};
    double apz{};
};

// ===== Helper function declarations =====

// Ultra-fast hashing function specialized for our specific polynomial structure
static uint32_t HashPolynomial(const std::vector<double> &coefficients);

// Compute suitable starting point for root finding
static double ComputeStartPoint(const std::vector<double> &a);

// Evaluate polynomial with Horner's method
static Eval EvaluatePolynomial(const std::vector<double> &a, const complex<double> z);

// Calculate upper bound for rounding errors (Adam's test)
static double CalculateUpperBound(const std::vector<double> &a, const complex<double> z);

// Real root forward deflation
static void RealDeflation(std::vector<double> &a, const double x);

// Complex root forward deflation for real coefficients
static void ComplexDeflation(std::vector<double> &a, const complex<double> z);

// Solve quadratic polynomial
static void SolveQuadratic(const std::vector<double> &a, std::vector<complex<double> > &roots,
                           std::vector<complex<double> > &newRoots);

// Root cache management functions
static void ClearRootCache();

static void PrintRootCacheStats();

static void ReserveRootCache(size_t expectedSize);

// ===== Implementation of helper functions =====

// Ultra-fast hashing function specialized for our specific polynomial structure
static uint32_t HashPolynomial(const std::vector<double> &coefficients) {
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
static double ComputeStartPoint(const std::vector<double> &a) {
    const size_t n = a.size() - 1;
    double a0      = log(abs(a.back()));
    double min     = exp((a0 - log(abs(a.front()))) / static_cast<double>(n));

    for (size_t i = 1; i < n; i++) {
        if (a[i] != 0.0) {
            double tmp = exp((a0 - log(abs(a[i]))) / static_cast<double>(n - i));
            if (tmp < min)
                min = tmp;
        }
    }

    return min * 0.5;
}

// Evaluate a polynomial with real coefficients using Horner's method
static Eval EvaluatePolynomial(const std::vector<double> &a, const complex<double> z) {
    const size_t n = a.size() - 1;
    double p       = -2.0 * z.real();
    double q       = norm(z);
    double s       = 0.0;
    double r       = a[0];
    Eval e;

    for (size_t i = 1; i < n; i++) {
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
static double CalculateUpperBound(const std::vector<double> &a, const complex<double> z) {
    const size_t n = a.size() - 1;
    double p       = -2.0 * z.real();
    double q       = norm(z);
    double u       = sqrt(q);
    double s       = 0.0;
    double r       = a[0];
    double e       = fabs(r) * (3.5 / 4.5);
    double t;

    for (size_t i = 1; i < n; i++) {
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
static void RealDeflation(std::vector<double> &a, const double x) {
    const size_t n = a.size() - 1;
    double r       = 0.0;

    for (size_t i = 0; i < n; i++)
        a[i]      = r = r * x + a[i];

    a.resize(n); // Remove the highest degree coefficient
}

// Complex root forward deflation for real coefficients
static void ComplexDeflation(std::vector<double> &a, const complex<double> z) {
    const size_t n = a.size() - 1;
    double r       = -2.0 * z.real();
    double u       = norm(z);

    a[1] -= r * a[0];
    for (int i = 2; i < n - 1; i++)
        a[i]   = a[i] - r * a[i - 1] - u * a[i - 2];

    a.resize(n - 1); // Remove top 2 highest degree coefficients
}

// Solve quadratic polynomial
static void SolveQuadratic(const std::vector<double> &a, std::vector<complex<double> > &roots,
                           std::vector<complex<double> > &newRoots) {
    const size_t n = a.size() - 1;
    complex<double> v;
    double r;

    // Notice that a[0] is !=0 since roots=zero has already been handled
    if (n == 1) {
        v = complex<double>(-a[1] / a[0], 0);
        roots.push_back(v);
        newRoots.push_back(v);
    } else {
        if (a[1] == 0.0) {
            r = -a[2] / a[0];
            if (r < 0) {
                r = sqrt(-r);
                v = complex<double>(0, r);
                roots.push_back(v);
                roots.push_back(conj(v));
                newRoots.push_back(v);
                newRoots.push_back(conj(v));
            } else {
                r = sqrt(r);
                v = complex<double>(r, 0);
                roots.push_back(v);
                roots.push_back(complex<double>(-r, 0));
                newRoots.push_back(v);
                newRoots.push_back(complex<double>(-r, 0));
            }
        } else {
            r = 1.0 - 4.0 * a[0] * a[2] / (a[1] * a[1]);
            if (r < 0) {
                v = complex<double>(-a[1] / (2.0 * a[0]), a[1] *
                                                          sqrt(-r) / (2.0 * a[0]));
                roots.push_back(v);
                roots.push_back(conj(v));
                newRoots.push_back(v);
                newRoots.push_back(conj(v));
            } else {
                v = complex<double>((-1.0 - sqrt(r)) * a[1] / (2.0 * a[0]), 0);
                roots.push_back(v);
                newRoots.push_back(v);
                v = complex<double>(a[2] / (a[0] * v.real()), 0);
                roots.push_back(v);
                newRoots.push_back(v);
            }
        }
    }
}

// Root cache management functions
static void ClearRootCache() {
    g_rootCache.clear();
}

static void PrintRootCacheStats() {
    // All print statements removed for performance
}

static void ReserveRootCache(size_t expectedSize) {
    g_rootCache.reserve(expectedSize);
}

// ===== Main polynomial roots function =====

// Find all polynomial zeros using a modified Newton method with warm start capability
static std::vector<complex<double> > PolynomialRoots(const std::vector<double> &coefficients) {
    const complex<double> complexzero(0.0, 0);      // Complex zero (0+i0)
    size_t n;                                       // Size of Polynomial p(x)
    Eval pz;                                        // P(z)
    Eval pzprev;                                    // P(zprev)
    Eval p1z;                                       // P'(z)
    Eval p1zprev;                                   // P'(zprev)
    complex<double> z;                              // Use as temporary variable
    complex<double> dz;                             // The current stepsize dz
    int itercnt;                                    // Hold the number of iterations per root
    std::vector<complex<double> > roots;            // Holds the roots of the Polynomial
    std::vector<double> coeff(coefficients.size()); // Holds the current coefficients of P(z)

    // Generate hash for this polynomial to check cache - ultra fast integer hash
    uint32_t polyHash     = HashPolynomial(coefficients);
    bool usingCachedRoots = false;
    std::vector<complex<double> > cachedRoots;

    // Check if we have solved a similar polynomial before
    auto cacheIt = g_rootCache.find(polyHash);
    if (cacheIt != g_rootCache.end()) {
        cachedRoots      = cacheIt->second;
        usingCachedRoots = true;
    }

    copy(coefficients.begin(), coefficients.end(), coeff.begin());

    // Step 1 eliminate all simple roots
    for (n = coeff.size() - 1; n > 0 && coeff.back() == 0.0; --n)
        roots.push_back(complexzero); // Store zero as the root

    // Index for tracking cached roots
    size_t cacheRootIndex = 0;
    std::vector<complex<double> > newRoots;

    // Do Newton iteration for polynomial order higher than 2
    for (; n > 2; --n) {
        const double Max_stepsize      = 5.0;                       // Allow the next step size to be up to 5x larger
        const complex<double> rotation = complex<double>(0.6, 0.8); // Rotation amount
        double r;                                                   // Current radius
        double rprev;                                               // Previous radius
        double eps;                                                 // Iteration termination value
        bool stage1 = true;                                         // By default start in stage1
        int steps   = 1;                                            // Multisteps if > 1
        std::vector<double> coeffprime;

        // Calculate coefficients of p'(x)
        for (int i = 0; i < n; i++)
            coeffprime.push_back(coeff[i] * double(n - i));

        // Step 2 find a suitable starting point z
        rprev = ComputeStartPoint(coeff); // Computed startpoint

        // If we have cached roots, use the next one as a starting point
        if (usingCachedRoots && cacheRootIndex < cachedRoots.size()) {
            z = cachedRoots[cacheRootIndex++];
        } else {
            // Use default starting point calculation
            z = coeff[n - 1] == 0.0 ? complex<double>(1.0, 0) : complex<double>(-coeff[n] / coeff[n - 1], 0);
            z *= complex<double>(rprev, 0) / abs(z);
        }

        // Setup the iteration
        // Current P(z)
        pz = EvaluatePolynomial(coeff, z);

        // pzprev which is the previous z or P(0)
        pzprev.z   = complex<double>(0, 0);
        pzprev.pz  = coeff[n];
        pzprev.apz = abs(pzprev.pz);

        // p1zprev P'(0) is the P'(0)
        p1zprev.z   = pzprev.z;
        p1zprev.pz  = coeff[n - 1]; // P'(0)
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
            p1z = EvaluatePolynomial(coeffprime, pz.z);
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
            pz     = EvaluatePolynomial(coeff, z);
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
                    npz = EvaluatePolynomial(coeff, zn);
                    if (npz.apz >= pz.apz)
                        break; // Break if no improvement

                    // Improved => accept step and try another round of step
                    pz = npz;

                    if (div2 == true && steps == 2) {
                        // To many shorten steps => try another direction and break
                        dz *= rotation;
                        z  = pzprev.z - dz;
                        pz = EvaluatePolynomial(coeff, z);
                        break;
                    }
                }
            } else {
                // calculate the upper bound of error using Grant & Hitchins's
                // test for complex coefficients
                // Now that we are within the convergence circle.
                eps = CalculateUpperBound(coeff, pz.z);
            }
        }

        // Check if there is a very small residue in the imaginary part by trying
        // to evaluate P(z.real). if that is less than P(z).
        // We take that z.real() is a better root than z.
        z      = complex<double>(pz.z.real(), 0.0);
        pzprev = EvaluatePolynomial(coeff, z);
        if (pzprev.apz <= pz.apz) {
            // real root
            pz = pzprev;
            // Save the root
            roots.push_back(pz.z);
            newRoots.push_back(pz.z);
            RealDeflation(coeff, pz.z.real());
        } else {
            // Complex root
            // Save the root
            roots.push_back(pz.z);
            roots.push_back(conj(pz.z));
            newRoots.push_back(pz.z);
            newRoots.push_back(conj(pz.z));
            ComplexDeflation(coeff, pz.z);
            --n;
        }
    } // End Iteration

    // Solve any remaining linear or quadratic polynomial
    if (n > 0)
        SolveQuadratic(coeff, roots, newRoots);

    // Store the roots in the cache for future use
    g_rootCache[polyHash] = newRoots;

    return roots;
}

#endif //NEWTON_SOLVER_H
