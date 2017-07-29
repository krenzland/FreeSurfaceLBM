#include "collision.hpp"
#include "freeSurface.hpp"

void computePostCollisionDistributions(double *currentCell, double tau, const double *const feq) {

    const double tau_1 = 1 / tau;

    for (int i = 0; i < Q; ++i) {
        currentCell[i] = currentCell[i] - tau_1 * (currentCell[i] - feq[i]);
    }
}

void doCollision(std::vector<double> &distributions, std::vector<double> &mass,
                 std::vector<double> &density, const std::vector<flag_t> &flagField, double tau,
                 const std::array<double, 3> &gravity, const coord_t &length, gridSet_t &filled,
                 gridSet_t &emptied) {
#pragma omp parallel for
    for (int z = 0; z < length[2] + 2; ++z) {
        double curDensity = 0;
        double velocity[3] = {0};
        double feq[19] = {0};

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

                // Set mass equal to density to avoid numerical instabilities.
                if (flagField[flagIndex] == flag_t::FLUID) {
                    mass[flagIndex] = curDensity;
                }
            }
        }
    }
}
