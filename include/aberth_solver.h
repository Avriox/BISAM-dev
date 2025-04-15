#ifndef ABERTH_SOLVER_H
#define ABERTH_SOLVER_H

/*
 *******************************************************************************
 *
 * Module name     :   AberthSolver.h
 * Module ID Nbr   :
 * Description     :   Solve n degree polynomial using Aberth-Ehrlich method
 * --------------------------------------------------------------------------
 * Change Record   :
 *
 * Version      Author/Date     Description of changes
 * -------     -------------  ----------------------
 * 01.01        ---/15Apr2025  Initial release
 * 01.02        ---/15Apr2025  Fixed warm start and convergence issues
 * 01.03        ---/15Apr2025  Ultra-optimized implementation using no complex numbers
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

#include <cmath>
#include <cstring>
#include <cstdlib>
#include "newton_solver.h" // Include Newton solver for shared definitions

using namespace std;

constexpr int MAX_ABERTH_ITER          = 100;
constexpr double OSCILLATION_THRESHOLD = 1e-10; // Threshold to detect oscillation
constexpr double GOOD_ENOUGH_THRESHOLD = 1e-5;  // Early termination threshold
constexpr int MAX_OSCILLATION_COUNT    = 4;     // Maximum oscillation count before early termination

// ===== Complex number structures and operations =====

// Structure representing a complex number without using std::complex
// struct Complex {
//     double re; // Real part
//     double im; // Imaginary part
//
//     // Fast norm (squared magnitude)
//     inline double norm() const {
//         return re * re + im * im;
//     }
//
//     // Fast absolute value (magnitude)
//     inline double abs() const {
//         return sqrt(norm());
//     }
// };
//
// // Structure for polynomial evaluation
// struct Eval {
//     double z_re;  // Real part of z
//     double z_im;  // Imaginary part of z
//     double pz_re; // Real part of P(z)
//     double pz_im; // Imaginary part of P(z)
//     double apz;   // |P(z)| - absolute value
// };
//
// // ===== Complex arithmetic operations =====
//
// // Complex number addition
// inline void complex_add(double &result_re, double &result_im,
//                         double a_re, double a_im,
//                         double b_re, double b_im) {
//     result_re = a_re + b_re;
//     result_im = a_im + b_im;
// }
//
// // Complex number subtraction
// inline void complex_sub(double &result_re, double &result_im,
//                         double a_re, double a_im,
//                         double b_re, double b_im) {
//     result_re = a_re - b_re;
//     result_im = a_im - b_im;
// }
//
// // Complex number multiplication
// inline void complex_mul(double &result_re, double &result_im,
//                         double a_re, double a_im,
//                         double b_re, double b_im) {
//     double tmp_re = a_re * b_re - a_im * b_im;
//     double tmp_im = a_re * b_im + a_im * b_re;
//     result_re     = tmp_re;
//     result_im     = tmp_im;
// }
//
// // Complex number division
// inline void complex_div(double &result_re, double &result_im,
//                         double a_re, double a_im,
//                         double b_re, double b_im) {
//     double denom  = b_re * b_re + b_im * b_im;
//     double tmp_re = (a_re * b_re + a_im * b_im) / denom;
//     double tmp_im = (a_im * b_re - a_re * b_im) / denom;
//     result_re     = tmp_re;
//     result_im     = tmp_im;
// }
//
// // Complex number absolute value (magnitude)
// inline double complex_abs(double re, double im) {
//     return sqrt(re * re + im * im);
// }
//
// // Complex number conjugate
// inline void complex_conj(double &result_re, double &result_im,
//                          double a_re, double a_im) {
//     result_re = a_re;
//     result_im = -a_im;
// }
//
// // Scale complex number by a real value
// inline void complex_scale(double &result_re, double &result_im,
//                           double a_re, double a_im,
//                           double scale) {
//     result_re = a_re * scale;
//     result_im = a_im * scale;
// }

// ===== Helper function declarations =====

// Helper function to eliminate zero roots (roots that are exactly zero)
static int ZeroRoots(int n, double a_re[], double a_im[], double res_re[], double res_im[]);

// Helper function to evaluate a polynomial at a complex point using Horner's method
static Eval Horner(int n, const double a_re[], const double a_im[], double z_re, double z_im);

// Helper function to calculate starting points for the algorithm
static void StartPoints(int n, const double apolyr[], double Z_re[], double Z_im[]);

// Helper function to solve quadratic and linear equations directly
static void Quadratic(int n, double a_re[], double a_im[], double res_re[], double res_im[]);

// ===== Implementation of helper functions =====

// Helper function to eliminate zero roots (roots that are exactly zero)
static int ZeroRoots(int n, double a_re[], double a_im[], double res_re[], double res_im[]) {
    int i = n;
    int k = 0;

    // Find roots that are exactly zero (coefficient a[n] = 0)
    while (i > 0 && complex_abs(a_re[i], a_im[i]) < 1e-14) {
        k++;
        i--;
    }

    // If we found zero roots, store them in the result array
    for (int j = n; j > i; j--) {
        res_re[j] = 0.0;
        res_im[j] = 0.0;
    }

    return i; // Return the new degree of the polynomial
}

// Helper function to evaluate a polynomial at a complex point using Horner's method
static Eval Horner(int n, const double a_re[], const double a_im[], double z_re, double z_im) {
    Eval result;
    double result_re = a_re[0];
    double result_im = a_im[0];
    double tmp_re, tmp_im;

    for (int i = 1; i <= n; i++) {
        // result = result * z + a[i]
        complex_mul(tmp_re, tmp_im, result_re, result_im, z_re, z_im);
        complex_add(result_re, result_im, tmp_re, tmp_im, a_re[i], a_im[i]);
    }

    result.z_re  = z_re;
    result.z_im  = z_im;
    result.pz_re = result_re;
    result.pz_im = result_im;
    result.apz   = complex_abs(result_re, result_im);

    return result;
}

// Helper function to calculate starting points for the algorithm
static void StartPoints(int n, const double apolyr[], double Z_re[], double Z_im[]) {
    // Calculate a suitable radius for the starting points
    double r0 = 0.0;
    for (int i = 0; i < n; i++) {
        if (apolyr[i] != 0.0) {
            double tmp = pow(apolyr[i] / apolyr[0], 1.0 / (n - i));
            if (tmp > r0) r0 = tmp;
        }
    }

    // Place starting points uniformly on a circle with radius r0
    const double PI = 3.14159265358979323846;
    r0 *= 1.1; // Slightly larger radius for better convergence

    for (int k = 1; k <= n; k++) {
        double angle = 2 * PI * (k - 1) / n;
        Z_re[k]      = r0 * cos(angle);
        Z_im[k]      = r0 * sin(angle);
    }
}

// Helper function to solve quadratic and linear equations directly
static void Quadratic(int n, double a_re[], double a_im[], double res_re[], double res_im[]) {
    if (n == 1) {
        // Linear equation: a[0]x + a[1] = 0
        complex_div(res_re[1], res_im[1], -a_re[1], -a_im[1], a_re[0], a_im[0]);
    } else if (n == 2) {
        // Quadratic equation: a[0]x^2 + a[1]x + a[2] = 0
        // discriminant = a[1]^2 - 4*a[0]*a[2]
        double disc_re, disc_im;
        double term1_re, term1_im;
        double term2_re, term2_im;
        double sqrt_disc_re, sqrt_disc_im;

        // Calculate a[1]^2
        complex_mul(term1_re, term1_im, a_re[1], a_im[1], a_re[1], a_im[1]);

        // Calculate 4*a[0]*a[2]
        complex_mul(term2_re, term2_im, a_re[0], a_im[0], a_re[2], a_im[2]);
        term2_re *= 4.0;
        term2_im *= 4.0;

        // Calculate discriminant = a[1]^2 - 4*a[0]*a[2]
        complex_sub(disc_re, disc_im, term1_re, term1_im, term2_re, term2_im);

        // Calculate sqrt(discriminant)
        // For complex square root, we use:
        // sqrt(r*e^(i*θ)) = sqrt(r)*e^(i*θ/2)
        double r     = complex_abs(disc_re, disc_im);
        double theta = atan2(disc_im, disc_re);
        sqrt_disc_re = sqrt(r) * cos(theta / 2);
        sqrt_disc_im = sqrt(r) * sin(theta / 2);

        // Calculate -a[1]
        double neg_a1_re = -a_re[1];
        double neg_a1_im = -a_im[1];

        // Calculate 2*a[0]
        double two_a0_re = 2.0 * a_re[0];
        double two_a0_im = 2.0 * a_im[0];

        // Calculate (-a[1] + sqrt(discriminant))/(2*a[0])
        double numer1_re, numer1_im;
        complex_add(numer1_re, numer1_im, neg_a1_re, neg_a1_im, sqrt_disc_re, sqrt_disc_im);
        complex_div(res_re[1], res_im[1], numer1_re, numer1_im, two_a0_re, two_a0_im);

        // Calculate (-a[1] - sqrt(discriminant))/(2*a[0])
        double numer2_re, numer2_im;
        complex_sub(numer2_re, numer2_im, neg_a1_re, neg_a1_im, sqrt_disc_re, sqrt_disc_im);
        complex_div(res_re[2], res_im[2], numer2_re, numer2_im, two_a0_re, two_a0_im);
    }
}

// ===== Main Aberth-Ehrlich function with anti-oscillation measures =====

// Main implementation of the Aberth-Ehrlich method
static void AberthEhrlich(int n, const double coeff_re[], const double coeff_im[],
                          double res_re[], double res_im[],
                          const double initialGuesses_re[] = nullptr,
                          const double initialGuesses_im[] = nullptr,
                          int numGuesses                   = 0) {
    bool dz_flag;
    int itercnt, i, j;
    double f, f0, f1, max_f, eps;

    // Complex variables
    double z_re, z_im, zi_re, zi_im, dz_re, dz_im;

    // Arrays for computation
    double *a_re, *a_im, *w_re, *w_im, *Z_re, *Z_im;
    bool *finish;
    double *apolyr;

    // Temporary variables for complex operations
    double tmp_re, tmp_im;

    // Evaluation results
    Eval fz, fz0, fz1;

    // Oscillation detection variables
    double prev_max_f       = 0.0;
    int oscillation_counter = 0;
    double damping_factor   = 1.0; // Start with no damping

    // Copy the original coefficients
    a_re = new double[n + 1];
    a_im = new double[n + 1];
    for (i = 0; i <= n; i++) {
        a_re[i] = coeff_re[i];
        a_im[i] = coeff_im[i];
    }

    // Eliminate zero roots
    n = ZeroRoots(n, a_re, a_im, res_re, res_im);

    if (n > 2) {
        double *a1_re = new double[n];
        double *a1_im = new double[n];

        // Calculate coefficients of f'(x)
        for (i = 0; i < n; i++) {
            a1_re[i] = a_re[i] * (n - i);
            a1_im[i] = a_im[i] * (n - i);
        }

        w_re   = new double[n + 1];
        w_im   = new double[n + 1];
        apolyr = new double[n + 1];
        Z_re   = new double[n + 1];
        Z_im   = new double[n + 1];
        finish = new bool[n + 1];

        // Simple upper bound for P(z) using Horner with complex coefficients
        f0  = complex_abs(a_re[n], a_im[n]);
        eps = 6 * n * f0 * pow((double) _DBL_RADIX, -DBL_MANT_DIG);

        for (i        = 0; i <= n; i++)
            apolyr[i] = complex_abs(a_re[i], a_im[i]);

        // If we have initial guesses (warm start), use them
        if (initialGuesses_re != nullptr && initialGuesses_im != nullptr && numGuesses >= n) {
            DEBUG_PRINT("Using warm start with " << numGuesses << " initial guesses");
            for (i = 1; i <= n; i++) {
                Z_re[i] = initialGuesses_re[i - 1]; // Adjust for 1-indexing in Z
                Z_im[i] = initialGuesses_im[i - 1];
            }
        } else {
            // Otherwise generate starting points
            StartPoints(n, apolyr, Z_re, Z_im);
        }

        for (i        = 1; i <= n; i++)
            finish[i] = false;

        max_f   = 1;
        dz_flag = true;

        // Start iteration
        for (itercnt = 1; dz_flag && max_f > eps && itercnt < MAX_ABERTH_ITER; itercnt++) {
            max_f   = 0;
            dz_flag = false;

            // Early termination for solutions that are "good enough"
            if (max_f < GOOD_ENOUGH_THRESHOLD && itercnt > 3) {
                DEBUG_PRINT("Early termination at iteration " << itercnt
                    << " - solution is good enough (error: " << max_f << ")");
                break;
            }

            // Escalated oscillation handling - if we've been oscillating too long, break out
            if (oscillation_counter > MAX_OSCILLATION_COUNT) {
                DEBUG_PRINT("Breaking out of persistent oscillation at iteration " << itercnt
                    << " - accepting current solution (error: " << max_f << ")");
                break;
            }

            // Check for oscillation
            if (itercnt > 2 && fabs(max_f - prev_max_f) < OSCILLATION_THRESHOLD) {
                oscillation_counter++;
                if (oscillation_counter > 3) {
                    // Apply stronger damping to break oscillation
                    damping_factor *= 0.5;
                    if (damping_factor < 0.01) damping_factor = 0.01; // Minimum damping

                    // Try random perturbation to break symmetry
                    if (oscillation_counter > 4) {
                        for (i = 1; i <= n; i++) {
                            if (!finish[i]) {
                                // Add a tiny random perturbation
                                double scale = 1e-15 * complex_abs(Z_re[i], Z_im[i]);
                                Z_re[i] += scale * ((double) rand() / RAND_MAX - 0.5);
                                Z_im[i] += scale * ((double) rand() / RAND_MAX - 0.5);
                            }
                        }
                        DEBUG_PRINT("Applied random perturbation to break oscillation");
                    } else {
                        DEBUG_PRINT("Oscillation detected, applying damping factor: " << damping_factor);
                    }
                }
            } else {
                // Reset oscillation counter if error is improving significantly
                if (prev_max_f > 0 && max_f < 0.5 * prev_max_f) {
                    oscillation_counter = 0;
                    // Gradually restore damping if not oscillating
                    if (damping_factor < 1.0) damping_factor *= 1.2;
                    if (damping_factor > 1.0) damping_factor = 1.0;
                }
            }

            prev_max_f = max_f;

            // Count how many roots are left to find
            int unfinished_roots = 0;
            for (i = 1; i <= n; i++) {
                if (!finish[i]) unfinished_roots++;
            }

            // If we're down to just a few roots, check if we can finish them directly
            if (unfinished_roots <= 2 && itercnt > 10) {
                // Try to solve directly for the remaining roots
                bool all_done = true;
                for (i = 1; i <= n; i++) {
                    if (!finish[i]) {
                        zi_re = Z_re[i];
                        zi_im = Z_im[i];
                        fz0   = Horner(n, a_re, a_im, zi_re, zi_im);
                        if (fz0.apz < 1e-10) {
                            finish[i] = true;
                            DEBUG_PRINT("Directly accepted root " << i << " as solution (error: " << fz0.apz << ")");
                        } else {
                            all_done = false;
                        }
                    }
                }
                if (all_done) break;
            }

            for (i = 1; i <= n; i++) {
                if (finish[i] == true) continue;

                zi_re = Z_re[i];
                zi_im = Z_im[i];
                fz0   = Horner(n, a_re, a_im, zi_re, zi_im);
                f0    = fz0.apz;
                fz1   = Horner(n - 1, a1_re, a1_im, zi_re, zi_im);
                f1    = fz1.apz;

                // Calculate w[i] = sum(1 / (zi - Z[j])) for j != i
                w_re[i] = 0.0;
                w_im[i] = 0.0;

                for (j = 1; j <= n; j++) {
                    if (i != j) {
                        // dz = zi - Z[j]
                        complex_sub(dz_re, dz_im, zi_re, zi_im, Z_re[j], Z_im[j]);

                        // tmp = 1 / dz
                        complex_div(tmp_re, tmp_im, 1.0, 0.0, dz_re, dz_im);

                        // w[i] += tmp
                        complex_add(w_re[i], w_im[i], w_re[i], w_im[i], tmp_re, tmp_im);
                    }
                }

                // dz = fz1 / fz0 - w[i]
                complex_div(tmp_re, tmp_im, fz1.pz_re, fz1.pz_im, fz0.pz_re, fz0.pz_im);
                complex_sub(dz_re, dz_im, tmp_re, tmp_im, w_re[i], w_im[i]);

                // dz = 1 / dz
                complex_div(dz_re, dz_im, 1.0, 0.0, dz_re, dz_im);

                // Apply damping factor to prevent oscillation
                dz_re *= damping_factor;
                dz_im *= damping_factor;

                // Store dz in w[i] for later use
                w_re[i] = dz_re;
                w_im[i] = dz_im;

                // z = zi - dz
                complex_sub(z_re, z_im, zi_re, zi_im, dz_re, dz_im);

                // Evaluate at new point
                fz = Horner(n, a_re, a_im, z_re, z_im);
                f  = fz.apz;

                // Add line search to ensure improvement
                if (f > f0) {
                    // If the step increases error, try more step sizes
                    double best_alpha = 1.0;
                    double best_f     = f;
                    double best_z_re  = z_re;
                    double best_z_im  = z_im;

                    // Try more step sizes (0.8, 0.6, 0.4, 0.2, 0.1, 0.01)
                    const double alphas[] = {0.8, 0.6, 0.4, 0.2, 0.1, 0.01};
                    for (const double alpha: alphas) {
                        // z_new = zi - alpha * dz
                        double z_new_re, z_new_im;
                        complex_scale(tmp_re, tmp_im, dz_re, dz_im, alpha);
                        complex_sub(z_new_re, z_new_im, zi_re, zi_im, tmp_re, tmp_im);

                        // Evaluate at new point
                        Eval fz_new  = Horner(n, a_re, a_im, z_new_re, z_new_im);
                        double f_new = fz_new.apz;

                        if (f_new < best_f) {
                            best_f     = f_new;
                            best_z_re  = z_new_re;
                            best_z_im  = z_new_im;
                            best_alpha = alpha;
                        }
                    }

                    // Use the best step size
                    z_re = best_z_re;
                    z_im = best_z_im;
                    f    = best_f;

                    // If we're taking a very small step, that's another sign of oscillation
                    if (best_alpha <= 0.1) {
                        oscillation_counter++;
                    }
                }

                Z_re[i] = z_re;
                Z_im[i] = z_im;

                // Check if dz is effectively zero (z + dz == z)
                dz_flag = dz_flag || (z_re + dz_re != z_re || z_im + dz_im != z_im);

                if (f > max_f)
                    max_f = f;

                if (f <= eps || (z_re + dz_re == z_re && z_im + dz_im == z_im)) {
                    double z0_re, z0_im;
                    finish[i] = true;

                    if (fabs(z_re) >= fabs(z_im)) {
                        z0_re = z_re;
                        z0_im = 0.0;
                    } else {
                        z0_re = 0.0;
                        z0_im = z_im;
                    }

                    Eval fz0 = Horner(n, a_re, a_im, z0_re, z0_im);

                    if (fz0.apz <= f) {
                        Z_re[i] = z_re = z0_re;
                        Z_im[i] = z_im = z0_im;
                    }
                }
            }

            DEBUG_PRINT("Iteration " << itercnt << " max error: " << max_f);
        }

        for (i = 1; i <= n; i++) {
            res_re[i] = Z_re[i];
            res_im[i] = Z_im[i];
        }

        delete[] finish;
        delete[] Z_re;
        delete[] Z_im;
        delete[] w_re;
        delete[] w_im;
        delete[] a1_re;
        delete[] a1_im;
        delete[] apolyr;
    } else {
        Quadratic(n, a_re, a_im, res_re, res_im);
    }

    delete[] a_re;
    delete[] a_im;
    return;
}

// ===== Interface functions to match the Newton solver =====

// C-style interface that directly writes to real and imaginary arrays
// Designed to match the interface of the Newton solver
static int PolynomialRootsAberth(
    const double *coefficient_vector_ptr,
    int degree,
    double *real_zero_vector_ptr,
    double *imaginary_zero_vector_ptr) {
    // Generate hash for this polynomial to check cache
    uint32_t polyHash     = HashPolynomial(coefficient_vector_ptr);
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

    // Convert real coefficients to separate real/imag arrays
    double complexCoeffs_re[MAX_COEFF];
    double complexCoeffs_im[MAX_COEFF];
    for (int i = 0; i <= degree; i++) {
        complexCoeffs_re[i] = coefficient_vector_ptr[i];
        complexCoeffs_im[i] = 0.0;
    }

    // Aberth implementation is 1-indexed
    double roots_re[MAX_ROOTS + 1];
    double roots_im[MAX_ROOTS + 1];

    // Call the Aberth-Ehrlich implementation, with warm start if available
    if (usingCachedRoots) {
        AberthEhrlich(degree, complexCoeffs_re, complexCoeffs_im, roots_re, roots_im,
                      cachedReal, cachedImag, cachedRootCount);
    } else {
        AberthEhrlich(degree, complexCoeffs_re, complexCoeffs_im, roots_re, roots_im);
    }

    // Convert results to the output arrays
    int rootCount = 0;
    double newReal[MAX_ROOTS];
    double newImag[MAX_ROOTS];

    for (int i = 1; i <= degree; i++) {
        real_zero_vector_ptr[rootCount]      = roots_re[i];
        imaginary_zero_vector_ptr[rootCount] = roots_im[i];

        newReal[rootCount] = roots_re[i];
        newImag[rootCount] = roots_im[i];

        rootCount++;
    }

    // Store the roots in the cache for future use
    RootStorage storage;
    memcpy(storage.real, newReal, sizeof(double) * rootCount);
    memcpy(storage.imag, newImag, sizeof(double) * rootCount);
    storage.count         = rootCount;
    g_rootCache[polyHash] = storage;

    return rootCount;
}

#endif //ABERTH_SOLVER_H
