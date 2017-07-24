#include "computeCellValues.hpp"
#include "LBDefinitions.hpp"
#include <assert.h>

double computeDensity(const double *const currentCell) {
    double density = 0.0;
#pragma omp simd
    for (int i = 0; i < Q; ++i) {
        density += currentCell[i];
    }
    return density;
}

void computeVelocity(const double *const currentCell, double density, double *velocity) {
    const int dimension = 3;

    // Zero array
    for (int i = 0; i < dimension; ++i) {
        velocity[i] = 0.0;
    }

#pragma omp simd
    for (int i = 0; i < Q; ++i) {
        for (int j = 0; j < dimension; ++j) {
            velocity[j] += currentCell[i] * LATTICEVELOCITIES[i][j];
        }
    }
// divide by density to get velocity
#pragma omp simd
    for (int k = 0; k < dimension; ++k) {
        velocity[k] /= density;
    }
}

void computeFeq(double density, const double *const velocity, double *feq) {
    const int dimension = 3;

    const double cs_pow2 = C_S * C_S;
    const double cs_2pow4 = 2 * cs_pow2 * cs_pow2;
    const double cs_2pow2 = 2 * cs_pow2;

    // First compute the norm of the velocity, as it is independent of i
    double u_dot_u = 0.0;

#pragma omp simd
    for (int i = 0; i < dimension; ++i) {
        u_dot_u += velocity[i] * velocity[i];
    }

    for (int i = 0; i < Q; ++i) {
        double ci_dot_u = 0.0;
#pragma omp simd
        for (int j = 0; j < dimension; ++j) {
            ci_dot_u += LATTICEVELOCITIES[i][j] * velocity[j];
        }
        const double ci_dot_u_pow2 = ci_dot_u * ci_dot_u;

        const double scaling = LATTICEWEIGHTS[i] * density;
        const double eq = 1.0 + ci_dot_u / cs_pow2 + ci_dot_u_pow2 / cs_2pow4 - u_dot_u / cs_2pow2;

        feq[i] = scaling * eq;
        assert(feq[i] >= 0.0 && feq[i] <= 1.1);
    }
}
