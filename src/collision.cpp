#include "collision.hpp"
#include "freeSurface.hpp"

void computePostCollisionDistributions(double *currentCell, double tau, const double *const feq) {

    const double tau_1 = 1 / tau;

    for (int i = 0; i < Q; ++i) {
        currentCell[i] = currentCell[i] - tau_1 * (currentCell[i] - feq[i]);
    }
}

void doCollision(std::vector<double> &distributions, const std::vector<double> &mass, const std::vector<flag_t> &flagField,
                 double tau, const coord_t &length, gridSet_t &filled, gridSet_t &emptied) {
    double density = 0;
    double velocity[3] = {0};
    double feq[19] = {0};

    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int x = 0; x < length[0] + 2; ++x) {
                const int flagIndex = indexForCell(x, y, z, length);
                if (flagField[flagIndex] != flag_t::FLUID || flagField[flagIndex] != flag_t::INTERFACE)
                    continue;

                const int distrIndex = flagIndex * Q;

                density = computeDensity(&distributions[distrIndex]);
                computeVelocity(&distributions[distrIndex], density, velocity);
                computeFeq(density, velocity, feq);
                computePostCollisionDistributions(&distributions[distrIndex], tau, feq);

                if (flagField[flagIndex] == flag_t::INTERFACE) {
                    // Check whether we have to convert the interface to an emptied or fluid cell.
                    // Doesn't actually update the flags but pushes them to a queue.
                    // We do this here so we do not have to calculate the density again.
                    const coord_t coord = {x, y, z};
                    getPotentialUpdates(coord, mass[flagIndex], density, filled, emptied);
                }
            }
        }
    }
}
