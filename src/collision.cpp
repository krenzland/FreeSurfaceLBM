#include "collision.hpp"
#include "freeSurface.hpp"

void computePostCollisionDistributions(double *currentCell, double tau, const double *const feq) {

    const double tau_1 = 1 / tau;

    for (int i = 0; i < Q; ++i) {
        currentCell[i] = currentCell[i] - tau_1 * (currentCell[i] - feq[i]);
    }
}

void doCollision(std::vector<double> &distributions, const std::vector<double> &mass,
                 std::vector<double> &density, const std::vector<flag_t> &flagField, double tau,
                 const coord_t &length, gridSet_t &filled, gridSet_t &emptied) {
    double curDensity = 0;
    double velocity[3] = {0};
    double feq[19] = {0};
    const double gravity[3] = {0, 0, -0.00015};

    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int x = 0; x < length[0] + 2; ++x) {
                const int flagIndex = indexForCell(x, y, z, length);
                if (flagField[flagIndex] != flag_t::FLUID &&
                    flagField[flagIndex] != flag_t::INTERFACE)
                    continue;

                const int distrIndex = flagIndex * Q;

                curDensity = computeDensity(&distributions[distrIndex]);
                density[flagIndex] = curDensity;
                computeVelocity(&distributions[distrIndex], curDensity, velocity);
                // apply gravity to velocity
                for (int i = 0; i < 3; ++i) {
                    velocity[i] += gravity[i] * tau;
                }
                computeFeq(curDensity, velocity, feq);
                computePostCollisionDistributions(&distributions[distrIndex], tau, feq);

                if (flagField[flagIndex] == flag_t::INTERFACE) {
                    // Check whether we have to convert the interface to an emptied or fluid cell.
                    // Doesn't actually update the flags but pushes them to a queue.
                    // We do this here so we do not have to calculate the density again.
                    const coord_t coord = {x, y, z};
                    getPotentialUpdates(coord, mass[flagIndex], curDensity, filled, emptied);
                }
            }
        }
    }
}
