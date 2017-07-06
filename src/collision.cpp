#include "collision.hpp"

void computePostCollisionDistributions(double *currentCell, const double *const tau,
                                       const double *const feq) {

    const double tau_1 = 1 / *tau;

    for (int i = 0; i < Q; ++i) {
        currentCell[i] = currentCell[i] - tau_1 * (currentCell[i] - feq[i]);
    }
}

void doCollision(double *collideField, flag_t *flagField, const double *const tau,
                 int number_cells) {
    double density = 0;
    double velocity[3] = {0};
    double feq[19] = {0};

    int offset = 0;

    for (int i = 0; i < number_cells; ++i) {
        if (flagField[i] != flag_t::FLUID)
            continue;

        offset = i * Q;

        density = computeDensity(collideField + offset);
        computeVelocity(collideField + offset, density, velocity);
        computeFeq(density, velocity, feq);
        computePostCollisionDistributions(collideField + offset, tau, feq);
    }
}

void doCollision(std::vector<double> &collideField, std::vector<flag_t> &flagField, double tau) {
    doCollision(collideField.data(), flagField.data(), &tau, (int)flagField.size());
}
