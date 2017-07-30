#include "streaming.hpp"
#include "LBMHelper.hpp"
#include "freeSurface.hpp"
#include <assert.h>
#include <iostream>

int neighbouring_fi_cell_index(int x, int y, int z, int fi, const coord_t &length) {
    int new_x = x - LATTICEVELOCITIES[fi][0];
    int new_y = y - LATTICEVELOCITIES[fi][1];
    int new_z = z - LATTICEVELOCITIES[fi][2];

    return indexForCell(new_x, new_y, new_z, length);
}

void doStreaming(const std::vector<double> &collideField, std::vector<double> &streamField,
                 const std::vector<double> &mass, std::vector<double> &density,
                 const coord_t &length, const std::vector<flag_t> &flagField,
                 std::vector<neighborhood_t> &neighborhood) {
#pragma omp parallel for schedule(static)
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
                    // 1. Interface cells have empty cells, with no valid distributions.
                    // 2. To preserve balance, we need to reconstruct distributions along the
                    // interface-normal.
                    const auto coord = coord_t{x, y, z};
                    // Density contains the densities of the previous timestep!
                    const auto normal =
                        computeSurfaceNormal(collideField, density, flagField, length, mass, coord);

                    bool hasFluidNeighbors = false;
                    bool hasEmptyNeighbors = false;

                    for (int i = 0; i < Q; ++i) {
                        const auto &vel = LATTICEVELOCITIES[i];
                        const int neighFlag =
                            indexForCell(coord_t{x - vel[0], y - vel[1], z - vel[2]}, length);
                        if (flagIndex == neighFlag)
                            continue;

                        const bool isEmptyAdjacent = flagField[neighFlag] == flag_t::EMPTY;

                        hasFluidNeighbors =
                            hasFluidNeighbors || (flagField[neighFlag] == flag_t::FLUID);
                        hasEmptyNeighbors =
                            hasEmptyNeighbors || flagField[neighFlag] == flag_t::EMPTY;

                        const int inv = inverseVelocityIndex(i);
                        const auto &invVelocity = LATTICEVELOCITIES[inv];
                        const double dotProduct = normal[0] * invVelocity[0] +
                                                  normal[1] * invVelocity[1] +
                                                  normal[2] * invVelocity[2];

                        const bool isNormalDirection = dotProduct > 0.0;

                        if (isEmptyAdjacent || isNormalDirection) {
                            // We need to reconstruct this distribution with eq. (4.5).
                            const double atmosphericPressure = 1.0;
                            // Note that we have to calculate the velocity of the time step before,
                            // hence the choice of distribution field.
                            const double curDensity = density[flagIndex];
                            std::array<double, 3> velocity;
                            computeVelocity(&collideField[fieldIndex], curDensity, velocity.data());
                            std::array<double, Q> feq;
                            computeFeq(atmosphericPressure, velocity.data(), feq.data());

                            streamField[fieldIndex + i] =
                                feq[inv] + feq[i] - collideField[fieldIndex + inv];
                        }
                    }
                    const bool isStandardCell = (hasEmptyNeighbors && hasFluidNeighbors);
                    if (isStandardCell) {
                        neighborhood[flagIndex] = neighborhood_t::STANDARD;
                    } else if (!hasFluidNeighbors) {
                        neighborhood[flagIndex] = neighborhood_t::NO_FLUID_NEIGHBORS;
                    } else if (!hasEmptyNeighbors) {
                        neighborhood[flagIndex] = neighborhood_t::NO_EMPTY_NEIGHBORS;
                    }
                }
            }
        }
    }
}