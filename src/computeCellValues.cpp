#include "computeCellValues.hpp"
#include "LBDefinitions.hpp"
#include "LBMHelper.hpp"
#include <assert.h>
#include <iostream>
#include <vector>

double computeDensity(const double *const currentCell) {
    double density = 0.0;
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

    for (int i = 0; i < Q; ++i) {
        for (int j = 0; j < dimension; ++j) {
            velocity[j] += currentCell[i] * LATTICEVELOCITIES[i][j];
        }
    }
    // divide by density to get velocity
    for (int k = 0; k < dimension; ++k) {
        velocity[k] /= density;
    }
}

double computeStressTensor(const std::vector<double> &distributions, double *feq, int cellIndex) {
    // We only calculate the norm of the tensor.
    double magnitude = 0.0;
    const int dims = 3;
    for (int alpha = 0; alpha < dims; ++alpha) {
        for (int beta = 0; beta < dims; ++beta) {
            // This is one entry of the 3 x 3 stress tensor.
            double elem = 0.0;
            for (int i = 0; i < Q; ++i) {
                const auto &vel = LATTICEVELOCITIES[i];
                elem += vel[alpha] * vel[beta] * (distributions[cellIndex + i] - feq[i]);
            }
            magnitude += elem * elem;
        }
    }
    return std::sqrt(magnitude);
}

double computeLocalRelaxationTime(double tau, double stressTensorNorm, double smagConstant) {
    const double viscosity = (tau - 0.5) / 3.0;
    const double smagSqr = smagConstant * smagConstant;
    const double stress =
        (std::sqrt(viscosity * viscosity + 18 * smagSqr * stressTensorNorm) - viscosity) /
        (6.0 * smagSqr);
    assert(stress >= 0.0); // Always increase viscosity!
    return 3 * (viscosity + smagSqr * stress) + 0.5;
}

void computeFeq(double density, const double *const velocity, double *feq) {
    const int dimension = 3;

    const double cs_pow2 = C_S * C_S;
    const double cs_2pow4 = 2 * cs_pow2 * cs_pow2;
    const double cs_2pow2 = 2 * cs_pow2;

    // First compute the norm of the velocity, as it is independent of i
    double u_dot_u = 0.0;

    for (int i = 0; i < dimension; ++i) {
        u_dot_u += velocity[i] * velocity[i];
    }

    for (int i = 0; i < Q; ++i) {
        double ci_dot_u = 0.0;

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
