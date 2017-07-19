#include "streaming.hpp"
#include "LBMHelper.hpp"
#include "freeSurface.hpp"
#include <assert.h>

int neighbouring_fi_cell_index(int x, int y, int z, int fi, const coord_t &length) {
    int new_x = x - LATTICEVELOCITIES[fi][0];
    int new_y = y - LATTICEVELOCITIES[fi][1];
    int new_z = z - LATTICEVELOCITIES[fi][2];

    return indexForCell(new_x, new_y, new_z, length);
}

void doStreaming(const std::vector<double> &collideField, std::vector<double> &streamField,
                 const std::vector<double> &mass, const std::vector<flag_t> &flagField,
                 const coord_t &length) {
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int x = 0; x < length[0] + 2; ++x) {
                const int flagIndex = indexForCell(x, y, z, length);
                const int fieldIndex = flagIndex * Q;
                if (flagField[flagIndex] == flag_t::FLUID ||
                    flagField[flagIndex] == flag_t::INTERFACE) {
                    // Standard streaming step.
                    for (int i = 0; i < Q; ++i) {
                        const int neighbour = neighbouring_fi_cell_index(x, y, z, i, length) * Q;
                        streamField[fieldIndex + i] = collideField[neighbour + i];

                        assert(streamField[fieldIndex + i] >= 0.0);
                    }
                }
                if (flagField[flagIndex] == flag_t::INTERFACE) {
                    // For interface cells we have to do some things differently.
                    // The second pass over the distributions makes things easier.
                    // We need to deal with the following things:
                    // 1. Interface cells have empty cells, with no valid distr ibutions.
                    // 2. To preserve balance, we need to reconstruct distributions along the
                    // interface-normal.
                    const auto coord = coord_t{x, y, z};
                    const auto normal = computeSurfaceNormal(collideField, mass, coord, length);

                    for (int i = 0; i < Q; ++i) {
                        const int neighbourFlag = neighbouring_fi_cell_index(x, y, z, i, length);
                        const int neighbourDistr = neighbourFlag * Q;

                        const int inv = inverseVelocityIndex(i);
                        const auto &invVelocity = LATTICEVELOCITIES[inv];
                        const double dotProduct = normal[0] * invVelocity[0] +
                                                  normal[1] * invVelocity[1] +
                                                  normal[2] * invVelocity[2];
                        const bool isNormalDirection = dotProduct > 0.0;

                        if (flagField[neighbourFlag] == flag_t::EMPTY || isNormalDirection) {
                            // We need to reconstruct this distribution with eq. (4.5).
                            const double atmosphericPressure = 1.0;
                            // Note that we have to calculate the velocity of the time step before,
                            // hence the choice
                            // of distribution field.
                            const double neighDensity =
                                computeDensity(&collideField[neighbourDistr]);
                            std::array<double, 3> neighVelocity;
                            computeVelocity(&collideField[neighbourDistr], neighDensity,
                                            neighVelocity.data());
                            std::array<double, 19> neighFeq;
                            computeFeq(atmosphericPressure, neighVelocity.data(), neighFeq.data());

                            // The paper uses a push-stream step, we use a pull-stream step.
                            // This is why we invert all fluid directions.
                            // TODO: Verify this claim and generally the correctness of the
                            // reconstruction.
                            streamField[fieldIndex + i] =
                                neighFeq[inv] + neighFeq[i] - collideField[inv];
                        }
                    }
                }
            }
        }
    }
}