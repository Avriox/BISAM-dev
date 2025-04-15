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
#include <iostream>
#include <functional>
#include <cfloat>
#include <map>
#include <string>
#include <sstream>
#include "PolynomialRootFinder.h"

//[BISAM] _DBL_RADIX was not defined
# define _DBL_RADIX FLT_RADIX

using namespace std;
constexpr int MAX_ITER = 50;

// Global cache for storing polynomial roots based on coefficient hash
static std::map<std::string, std::vector<std::complex<double> > > g_rootCache;

// Helper function to generate a hash string from polynomial coefficients
// The precision parameter controls how many significant digits to use for the hash
static std::string HashPolynomial(const std::vector<double> &coefficients, int precision = 3) {
    std::ostringstream hashStream;
    hashStream.precision(precision);

    // We'll include the polynomial degree in the hash
    hashStream << coefficients.size() - 1 << ":";

    // Add the coefficients with limited precision to the hash
    for (const double &coef: coefficients) {
        hashStream << coef << ",";
    }

    return hashStream.str();
}

// Find all polynomial zeros using a modified Newton method with warm start capability
// 1) Check if we have previously solved a similar polynomial and use those roots
// 2) Eliminate all simple roots (roots equal to zero)
// 3) Find a suitable starting point or use cached roots
// 4) Find a root using Newton
// 5) Divide the root up in the polynomial reducing its order with one
// 6) Repeat steps 3 to 5 until the polynomial is of the order of two
//    whereafter the remaining one or two roots are found by the direct formula
// 7) Cache the found roots for future use
// Notice:
//      The coefficients for p(x) is stored in descending order. coefficients[0]
//      is a(n)x^n, coefficients[1] is a(n-1)x^(n-1),...,
//      coefficients[n-1] is a(1)x, coefficients[n] is a(0)
//
static std::vector<complex<double> > PolynomialRoots(const std::vector<double> &coefficients) {
    struct eval {
        complex<double> z{};
        complex<double> pz{};
        double apz{};
    };
    const complex<double> complexzero(0.0, 0);      // Complex zero (0+i0)
    size_t n;                                       // Size of Polynomial p(x)
    eval pz;                                        // P(z)
    eval pzprev;                                    // P(zprev)
    eval p1z;                                       // P'(z)
    eval p1zprev;                                   // P'(zprev)
    complex<double> z;                              // Use as temporary variable
    complex<double> dz;                             // The current stepsize dz
    int itercnt;                                    // Hold the number of iterations per root
    std::vector<complex<double> > roots;            // Holds the roots of the Polynomial
    std::vector<double> coeff(coefficients.size()); // Holds the current coefficients of P(z)

    // Generate hash for this polynomial to check cache
    std::string polyHash  = HashPolynomial(coefficients, 1);
    bool usingCachedRoots = false;
    std::vector<complex<double> > cachedRoots;

    // Check if we have solved a similar polynomial before
    auto cacheIt = g_rootCache.find(polyHash);
    if (cacheIt != g_rootCache.end()) {
        cachedRoots      = cacheIt->second;
        usingCachedRoots = true;
        // std::cout << "Found cached roots for polynomial with hash: " << polyHash << std::endl;
    }

    copy(coefficients.begin(), coefficients.end(), coeff.begin());
    // Step 1 eliminate all simple roots
    for (n = coeff.size() - 1; n > 0 && coeff.back() == 0.0; --n)
        roots.push_back(complexzero); // Store zero as the root

    // Compute the next starting point based on the polynomial coefficients
    // A root will always be outside the circle from the origin and radius min
    auto startpoint = [&](const std::vector<double> &a) {
        const size_t n = a.size() - 1;
        double a0      = log(abs(a.back()));
        double min     = exp((a0 - log(abs(a.front()))) / static_cast<double>(n));

        for (size_t i = 1; i < n; i++)
            if (a[i] != 0.0) {
                double tmp = exp((a0 - log(abs(a[i]))) / static_cast<double>(n - i));
                if (tmp < min)
                    min = tmp;
            }

        return min * 0.5;
    };

    // Evaluate a polynomial with real coefficients a[] at a complex point z and
    // return the result
    // This is Horner's method of avoiding complex arithmetic
    auto horner = [](const std::vector<double> &a, const complex<double> z) {
        const size_t n = a.size() - 1;
        double p       = -2.0 * z.real();
        double q       = norm(z);
        double s       = 0.0;
        double r       = a[0];
        eval e;

        for (size_t i = 1; i < n; i++) {
            double t = a[i] - p * r - q * s;
            s        = r;
            r        = t;
        }

        e.z   = z;
        e.pz  = complex<double>(a[n] + z.real() * r - q * s, z.imag() * r);
        e.apz = abs(e.pz);
        return e;
    };

    // Calculate a upper bound for the rounding errors performed in a
    // polynomial with real coefficient a[] at a complex point z.
    // (Adam's test)
    auto upperbound = [](const std::vector<double> &a, const complex<double> z) {
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
    };

    // Real root forward deflation.
    auto realdeflation = [&](std::vector<double> &a, const double x) {
        const size_t n = a.size() - 1;
        double r       = 0.0;

        for (size_t i = 0; i < n; i++)
            a[i]      = r = r * x + a[i];

        a.resize(n); // Remove the highest degree coefficients.
    };

    // Complex root forward deflation for real coefficients
    auto complexdeflation = [&](std::vector<double> &a, const complex<double> z) {
        const size_t n = a.size() - 1;
        double r       = -2.0 * z.real();
        double u       = norm(z);

        a[1] -= r * a[0];
        for (int i = 2; i < n - 1; i++)
            a[i]   = a[i] - r * a[i - 1] - u * a[i - 2];

        a.resize(n - 1); // Remove top 2 highest degree coefficienst
    };

    // Index for tracking cached roots
    size_t cacheRootIndex = 0;
    std::vector<complex<double> > newRoots;

    // Do Newton iteration for polynomial order higher than 2
    for (; n > 2; --n) {
        const double Max_stepsize = 5.0; // Allow the next step size to be
        // up to 5 times larger than the previous step size
        const complex<double> rotation = complex<double>(0.6, 0.8); // Rotation amount
        double r;                                                   // Current radius
        double rprev;                                               // Previous radius
        double eps;                                                 // The iteration termination value
        bool stage1 = true;                                         // By default it start the iteration in stage1
        int steps   = 1;                                            // Multisteps if > 1
        std::vector<double> coeffprime;

        // Calculate coefficients of p'(x)
        for (int i = 0; i < n; i++)
            coeffprime.push_back(coeff[i] * double(n - i));

        // Step 2 find a suitable starting point z
        rprev = startpoint(coeff); // Computed startpoint

        // If we have cached roots, use the next one as a starting point
        if (usingCachedRoots && cacheRootIndex < cachedRoots.size()) {
            z = cachedRoots[cacheRootIndex++];
            // std::cout << "Using cached root as starting point: " << z << std::endl;
        } else {
            // Use default starting point calculation
            z = coeff[n - 1] == 0.0 ? complex<double>(1.0, 0) : complex<double>(-coeff[n] / coeff[n - 1], 0);
            z *= complex<double>(rprev, 0) / abs(z);
        }

        // Setup the iteration
        // Current P(z)
        pz = horner(coeff, z);

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
        r  = rprev *= Max_stepsize; // Make a reasonable radius of the
        // maximum step size allowed
        // Preliminary eps computed at point P(0) by a crude estimation
        eps = 2 * n * pzprev.apz * pow((double) _DBL_RADIX, -DBL_MANT_DIG);

        // Start iteration and stop if z doesnt change or apz <= eps
        // we do z+dz!=z instead of dz!=0. if dz does not change z
        // then we accept z as a root
        for (itercnt = 0; pz.z + dz != pz.z &&
                          pz.apz > eps && itercnt < MAX_ITER; itercnt++) {
            // Calculate current P'(z)
            p1z = horner(coeffprime, pz.z);
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
            z      = pzprev.z - dz;    // Next z
            pz     = horner(coeff, z); //ff = pz.apz;
            steps  = 1;
            if (stage1) {
                // Try multiple steps or shorten steps depending if P(z)
                // is an improvement or not P(z)<P(zprev)
                bool div2;
                complex<double> zn;
                eval npz;

                zn = pz.z;
                for (div2 = pz.apz > pzprev.apz; steps <= n; ++steps) {
                    if (div2 == true) {
                        // Shorten steps
                        dz *= complex<double>(0.5, 0);
                        zn = pzprev.z - dz;
                    } else
                        zn -= dz; // try another step in the same direction

                    // Evaluate new try step
                    npz = horner(coeff, zn);
                    if (npz.apz >= pz.apz)
                        break; // Break if no improvement

                    // Improved => accept step and try another round of step
                    pz = npz;

                    if (div2 == true && steps == 2) {
                        // To many shorten steps => try another direction and break
                        dz *= rotation;
                        z  = pzprev.z - dz;
                        pz = horner(coeff, z);
                        break;
                    }
                }
            } else {
                // calculate the upper bound of error using Grant & Hitchins's
                // test for complex coefficients
                // Now that we are within the convergence circle.
                eps = upperbound(coeff, pz.z);
            }
        }

        // Print the number of iterations it took to converge for this root
        // std::cout << "Root found after " << itercnt << " iterations." << std::endl;

        // Check if warm start improved convergence
        // if (usingCachedRoots && cacheRootIndex <= cachedRoots.size()) {
        //     std::cout << "Using warm start: root found in " << itercnt
        //             << " iterations (vs. typically more without warm start)" << std::endl;
        // }

        // Check if there is a very small residue in the imaginary part by trying
        // to evaluate P(z.real). if that is less than P(z).
        // We take that z.real() is a better root than z.
        z      = complex<double>(pz.z.real(), 0.0);
        pzprev = horner(coeff, z);
        if (pzprev.apz <= pz.apz) {
            // real root
            pz = pzprev;
            // Save the root
            roots.push_back(pz.z);
            newRoots.push_back(pz.z);
            realdeflation(coeff, pz.z.real());
        } else {
            // Complex root
            // Save the root
            roots.push_back(pz.z);
            roots.push_back(conj(pz.z));
            newRoots.push_back(pz.z);
            newRoots.push_back(conj(pz.z));
            complexdeflation(coeff, pz.z);
            --n;
        }
    } // End Iteration

    // Solve any remaining linear or quadratic polynomial
    // For Polynomial with real coefficients a[],
    // The complex solutions are stored in the back of the roots
    auto quadratic = [&](const std::vector<double> &a) {
        const size_t n = a.size() - 1;
        complex<double> v;
        double r;

        // Notice that a[0] is !=0 since roots=zero has already been handle
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
        return;
    };

    if (n > 0)
        quadratic(coeff);

    // Store the roots in the cache for future use
    g_rootCache[polyHash] = newRoots;

    return roots;
}

// Helper function to clear the root cache if needed
static void ClearRootCache() {
    g_rootCache.clear();
    std::cout << "Root cache cleared" << std::endl;
}

// Helper function to get cache stats
static void PrintRootCacheStats() {
    std::cout << "Root cache contains " << g_rootCache.size() << " entries" << std::endl;
}

#endif //NEWTON_SOLVER_H
