#include "freeSurface.hpp"
#include <iostream>
#include <numeric>

double computeFluidFraction(const std::vector<double> &fluidFraction, const std::vector<flag_t> &flags, int idx,
                            const std::vector<double> &mass) {
    if (flags[idx] == flag_t::EMPTY) {
        return 0.0;
    } else if (flags[idx] == flag_t::INTERFACE) {
        const double vof = fluidFraction[idx];
        // in some cases it can be negative or larger than 1, e.g. if the cell is recently filled.
        return std::min(1.0, std::max(0.0, vof));
    } else {
        return 1.0;
    }
}

std::array<double, 3> computeSurfaceNormal(const std::vector<double> &fluidFraction, const coord_t &length, const std::vector<double> &mass,
                                           const coord_t &position, const std::vector<flag_t> &flags) {
    auto normal = std::array<double, 3>();
    // We approximate the surface normal element-wise using central-differences of the fluid
    // fraction gradient.
    for (size_t dim = 0; dim < normal.size(); ++dim) {
        coord_t curPosition = position;

        // We need the next and previous neighbour for dimension dim.
        curPosition[dim]++;
        const auto nextNeighbour = indexForCell(curPosition, length);
        curPosition[dim] -= 2;
        const auto prevNeighbour = indexForCell(curPosition, length);

        const double plusFluidFraction = computeFluidFraction(fluidFraction, flags, nextNeighbour, mass);
        const double minusFluidFraction = computeFluidFraction(fluidFraction, flags, prevNeighbour, mass);

        normal[dim] = 0.5 * (minusFluidFraction - plusFluidFraction);
    }
    return normal;
}

void streamMass(const std::vector<double> &distributions, std::vector<double> &fluidFraction, const coord_t &length,
                std::vector<double> &mass, const std::vector<neighborhood_t> &neighborhood,
                const std::vector<flag_t> &flags) {
#pragma omp parallel for schedule(static)
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int x = 0; x < length[0] + 2; ++x) {
                double deltaMass = 0.0;

                const coord_t curCell = coord_t{x, y, z};
                const int flagIndex = indexForCell(curCell, length);
                // We only consider the mass going into interface cells.
                // Empty cells have zero mass, full cells have mass 1.
                if (flags[flagIndex] != flag_t::INTERFACE)
                    continue;

                const int fieldIndex = flagIndex * Q;
                const double curFluidFraction =
                        computeFluidFraction(fluidFraction, flags, flagIndex, mass);

                for (int i = 0; i < Q; ++i) {
                    const auto &vel = LATTICEVELOCITIES[i];
                    const auto neighCell = coord_t{x + vel[0], y + vel[1], z + vel[2]};
                    const auto neighFlag = indexForCell(neighCell, length);
                    const auto neighField = neighFlag * Q;
                    if (neighFlag == flagIndex)
                        continue;

                    if (flags[neighFlag] == flag_t::FLUID) {
                        // Exchange interface and fluid at x + \Delta t e_i (eq. 4.2)
                        deltaMass += distributions[neighField + inverseVelocityIndex(i)] -
                                     distributions[fieldIndex + i];
                    } else if (flags[neighFlag] == flag_t::INTERFACE) {
                        const double neighFluidFraction =
                                computeFluidFraction(fluidFraction, flags, neighFlag, mass);
                        // Exchange interface and interface at x + \Delta t e_i (eq. 4.2)
                        const double s_e =
                            calculateSE(distributions, flags, curCell, length, i, neighborhood);
                        deltaMass += s_e * 0.5 * (curFluidFraction + neighFluidFraction);
                    }
                }
                mass[flagIndex] += deltaMass; // (eq. 4.4)
            }
        }
    }
}

double calculateSE(const std::vector<double> &distributions, const std::vector<flag_t> &flags,
                   const coord_t &coord, const coord_t &length, const int curFiIndex,
                   const std::vector<neighborhood_t> &neighborhood) {
    // Returns a slightly different mass exchange depending on the neighbourhood of the cells.
    // This reduces the amount of abandoned interface cells.

    const auto &vel = LATTICEVELOCITIES[curFiIndex];
    const auto neigh = coord_t{coord[0] + vel[0], coord[1] + vel[1], coord[2] + vel[2]};

    const int curFlag = indexForCell(coord, length);
    const int neighFlag = indexForCell(neigh, length);

    const int curCell = curFlag * Q;
    const int neighCell = neighFlag * Q;

    const int inv = inverseVelocityIndex(curFiIndex);

    if (neighborhood[curFlag] == neighborhood[neighFlag]) {
        return distributions[neighCell + inv] - distributions[curCell + curFiIndex];
    }

    if ((neighborhood[neighFlag] == neighborhood_t::STANDARD &&
         neighborhood[curFlag] == neighborhood_t::NO_FLUID_NEIGHBORS) ||
        (neighborhood[neighFlag] == neighborhood_t::NO_EMPTY_NEIGHBORS &&
         (neighborhood[curFlag] == neighborhood_t::STANDARD ||
          neighborhood[curFlag] == neighborhood_t::NO_FLUID_NEIGHBORS))) {
        return -distributions[neighCell + curFiIndex];
    }

    // Otherwise: (neigh == STANDARD && cur = NO_EMPTY) || (neigh = NO_FLUID && (cur = STANDARD ||
    // cur NO_EMPTY)
    // TODO FIX!
    return distributions[neighCell + inv];
}

void getPotentialUpdates(const std::vector<double> &mass, std::vector<double> &fluidFraction, std::vector<flag_t> &flags,
                         std::vector<neighborhood_t> &neighborhood, const coord_t &length) {
    // Check whether we have to convert the interface to an emptied or fluid cell.
    // Doesn't actually update the flags but pushes them to a queue.
    // We do this here so we do not have to calculate the density again.

    // Offset avoids periodically switching between filled and empty status.
    const double offset = 10e-3;
#pragma omp parallel for schedule(guided)
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int x = 0; x < length[0] + 2; ++x) {
                auto coord = coord_t{x, y, z};
                const int flagIndex = indexForCell(coord, length);
                if (flags[flagIndex] != flag_t::INTERFACE)
                    continue;
                // Eq. 4.7
                double curDensity;
                if (fluidFraction[flagIndex] != 0.0) {
                    curDensity = mass[flagIndex] / fluidFraction[flagIndex];
                } else {
                    curDensity = 0.0;
                }
                if (mass[flagIndex] > (1 + offset) * curDensity) {
                    flags[flagIndex] = flag_t::INTERFACE_TO_FLUID;
                } else if (mass[flagIndex] < -offset * curDensity) {
                    // Emptied
                    flags[flagIndex] = flag_t::INTERFACE_TO_EMPTY;
                }
                if (neighborhood[flagIndex] == neighborhood_t::NO_FLUID_NEIGHBORS &&
                    mass[flagIndex] < (0.1 * curDensity)) {
                    flags[flagIndex] = flag_t::INTERFACE_TO_EMPTY;
                }
                if (neighborhood[flagIndex] == neighborhood_t::NO_EMPTY_NEIGHBORS &&
                    mass[flagIndex] > (0.9 * curDensity)) {
                    flags[flagIndex] = flag_t::INTERFACE_TO_FLUID;
                }
            }
        }
    }
}

void interpolateEmptyCell(std::vector<double> &distributions, std::vector<double> &mass, std::vector<double> &fluidFraction,
                          const coord_t &coord, const coord_t &length, const std::vector<flag_t> &flags) {
    // Note: We only interpolate cells that are not emptied cells themselves!
    const int flagIndex = indexForCell(coord, length);
    const int cellIndex = flagIndex * Q;

    int numNeighs = 0;
    double avgDensity = 0.0;
    auto avgVel = std::array<double, 3>{};
    for (int i = 0; i < Q; ++i) {
        const auto &vel = LATTICEVELOCITIES[i];
        const auto neigh = coord_t{coord[0] + vel[0], coord[1] + vel[1], coord[2] + vel[2]};

        const int neighFlagIndex = indexForCell(neigh, length);
        if (neighFlagIndex == flagIndex)
            continue; // Ignore current cell!

        // Don't interpolate cells that are emptied or filled.
        if (flags[neighFlagIndex] == flag_t::FLUID || flags[neighFlagIndex] == flag_t::INTERFACE ||
            flags[neighFlagIndex] == flag_t::INTERFACE_TO_FLUID) {
            const int neighDistrIndex = neighFlagIndex * Q;
            const double neighDensity = computeDensity(&distributions[neighDistrIndex]);
            std::array<double, 3> neighVelocity;
            computeVelocity(&distributions[neighDistrIndex], neighDensity, neighVelocity.data());

            ++numNeighs;
            avgDensity += neighDensity;
            avgVel[0] += neighVelocity[0];
            avgVel[1] += neighVelocity[1];
            avgVel[2] += neighVelocity[2];
        }
    }

    // Every former empty cell has at least one interface cell as neighbour, otherwise we have a
    // worse problem than division by zero.
    assert(numNeighs != 0);
    avgDensity /= numNeighs;
    avgVel[0] /= numNeighs;
    avgVel[1] /= numNeighs;
    avgVel[2] /= numNeighs;

    fluidFraction[flagIndex] = mass[flagIndex] / avgDensity;
    auto feq = std::array<double, 19>{};
    computeFeq(avgDensity, avgVel.data(), feq.data());
    for (int i = 0; i < Q; ++i) {
        distributions[cellIndex + i] = feq[i];
    }
}

void flagReinit(std::vector<double> &distributions, std::vector<double> &mass, std::vector<double> &fluidFraction,
                const coord_t &length, std::vector<flag_t> &flags) {
// We set up the interface cells around all converted cells.
// It's very important to do this in a way such that the order of conversion doesn't matter.

// First set interface for all filled cells.
#pragma omp parallel for schedule(guided)
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = 1; x < length[0] + 1; ++x) {
                const int curFlag = indexForCell(x, y, z, length);
                if (flags[curFlag] != flag_t::INTERFACE_TO_FLUID)
                    continue;
                // Find all neighbours of this cell.
                for (const auto &vel : LATTICEVELOCITIES) {
                    coord_t neighbor = {x, y, z};
                    neighbor[0] += vel[0];
                    neighbor[1] += vel[1];
                    neighbor[2] += vel[2];
                    const int neighFlag = indexForCell(neighbor, length);
                    if (curFlag == neighFlag)
                        continue;

                    // This neighbor is converted to an interface cell iff. it is an empty cell or a
                    // cell that would
                    // become an emptied cell.
                    // We need to remove it from the emptied set, otherwise we might have holes in
                    // the interface.
                    if (flags[neighFlag] == flag_t::EMPTY) {
                        flags[neighFlag] = flag_t::INTERFACE;
                        mass[neighFlag] = 0.0;
                        fluidFraction[neighFlag] = 0.0;
                        // Notice that the new interface cells don't have any valid distributions.
                        // They are initialised with f^{eq}_i (p_{avg}, v_{avg}), which are the
                        // average density and
                        // velocity of all neighbouring fluid and interface cells.
                        interpolateEmptyCell(distributions, mass, fluidFraction, neighbor, length,
                                             flags);
                    } else if (flags[neighFlag] == flag_t::INTERFACE_TO_EMPTY) {
                        // Already is an interface but should not be converted to an empty cell
                        // later.
                        flags[neighFlag] = flag_t::INTERFACE;
                    }
                }
            }
        }
    }

// Now we can consider all filled cells!
#pragma omp parallel for schedule(guided)
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = 1; x < length[0] + 1; ++x) {
                const int curFlag = indexForCell(x, y, z, length);
                // Find all neighbours of this cell.
                if (flags[curFlag] != flag_t::INTERFACE_TO_EMPTY)
                    continue;
                for (const auto &vel : LATTICEVELOCITIES) {
                    coord_t neighbor = {x, y, z};
                    neighbor[0] += vel[0];
                    neighbor[1] += vel[1];
                    neighbor[2] += vel[2];
                    const int neighFlag = indexForCell(neighbor, length);
                    if (curFlag == neighFlag)
                        continue;

                    // This neighbor is converted to an interface cell iff. it is an empty cell or a
                    // cell that would
                    // become an emptied cell.
                    // We need to remove it from the emptied set, otherwise we might have holes in
                    // the interface.
                    if (flags[neighFlag] == flag_t::FLUID) {
                        flags[neighFlag] = flag_t::INTERFACE;
                        mass[neighFlag] = computeDensity(&distributions[neighFlag * Q]);
                        fluidFraction[neighFlag] = 1.0;
                        // We can reuse the distributions as they are still valid.
                    }
                }
            }
        }
    }
}

enum class update_t { FILLED, EMPTIED };

void distributeSingleMass(const std::vector<double> &distributions, std::vector<double> &mass,
                          std::vector<double> &massChange, const coord_t &length, const coord_t &coord,
                          std::vector<flag_t> &flags, std::vector<double> &fluidFraction, const update_t &type) {
    // First determine how much mass needs to be redistributed and fix mass of converted cell.
    const int flagIndex = indexForCell(coord, length);
    double excessMass;
    const double density = computeDensity(&distributions[flagIndex * Q]);
    if (type == update_t::FILLED) {
        // Interface -> Full cell, filled cells have mass and should have a mass equal to their
        // density.
        excessMass = mass[flagIndex] - density;

        #pragma omp atomic
        massChange[flagIndex] -= excessMass;
    } else {
        // Interface -> Empty cell, empty cells should not have any mass so all mass is excess mass.
        excessMass = mass[flagIndex];
        #pragma omp atomic
        massChange[flagIndex] -= excessMass;
    }

    /* The distribution of excess mass is surprisingly non-trivial.
       For a more detailed description, refer to pages 32f. of Thürey's thesis but here's the gist:
       We do not distribute the mass uniformly to all neighbouring interface cells but rather
       correct for balance.
       The reason for this is that the fluid interface moved beyond the current cell.
       We rebalance things by weighting the mass updates according to the direction of the interface
       normal.
       This has to be done in two steps, we first calculate all updated weights, normalize them and
       then, in a second step, update the weights.*/

    // Step 1: Calculate the unnormalized weights.
    const auto normal = computeSurfaceNormal(fluidFraction, length, mass, coord, flags);
    std::array<double, 19> weights{};
    std::array<double, 19> weightsBackup{}; // Sometimes first weights is all zero.

    for (size_t i = 0; i < LATTICEVELOCITIES.size(); ++i) {
        const auto &vel = LATTICEVELOCITIES[i];
        coord_t neighbor = {coord[0] + vel[0], coord[1] + vel[1], coord[2] + vel[2]};

        const int neighFlag = indexForCell(neighbor, length);
        if (flagIndex == neighFlag || flags[neighFlag] != flag_t::INTERFACE)
            continue;

        weightsBackup[i] = 1.0;

        const double dotProduct = normal[0] * vel[0] + normal[1] * vel[1] + normal[2] * vel[2];
        if (type == update_t::FILLED) {
            weights[i] = std::max(0.0, dotProduct);
        } else { // EMPTIED
            weights[i] = -std::min(0.0, dotProduct);
        }
        assert(weights[i] >= 0.0);
    }

    // Step 2: Calculate normalizer (otherwise sum of weights != 1.0)
    double normalizer = std::accumulate(weights.begin(), weights.end(), 0.0, std::plus<double>());

    if (normalizer == 0.0) {
        // Mass cannot be distributed along the normal, distribute to all interface cells equally.
        weights = weightsBackup;
        normalizer = std::accumulate(weights.begin(), weights.end(), 0.0, std::plus<double>());
    }
    if (normalizer == 0.0) {
        return;
    }

    // Step 3: Redistribute weights. As non-interface cells have weight 0, we can just loop through
    // all cells.
    for (size_t i = 0; i < LATTICEVELOCITIES.size(); ++i) {
        const auto &vel = LATTICEVELOCITIES[i];
        coord_t neighbor = {coord[0] + vel[0], coord[1] + vel[1], coord[2] + vel[2]};

        const int neighFlag = indexForCell(neighbor, length);

        #pragma omp atomic
        massChange[neighFlag] += (weights[i] / normalizer) * excessMass;
    }
}

void distributeMass(const std::vector<double> &distributions, std::vector<double> &mass, const coord_t &length,
                    std::vector<flag_t> &flags, std::vector<double> &fluidFraction) {
    // Here we redistribute the excess mass of the cells.
    // It is important that we get a copy of the filled/emptied where all converted cells are stored
    // and no other cells.
    // This excludes emptied cells that are used as interface cells instead!
    auto massChange = std::vector<double>(mass.size(), 0.0);

#pragma omp parallel for schedule(guided)
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = 1; x < length[0] + 1; ++x) {
                const int curFlag = indexForCell(x, y, z, length);
                const coord_t coord = {x, y, z};
                if (flags[curFlag] == flag_t::INTERFACE_TO_FLUID) {
                    distributeSingleMass(distributions, mass, massChange, length, coord,
                                         flags, fluidFraction, update_t::FILLED);
                    flags[curFlag] = flag_t::FLUID;
                } else if (flags[curFlag] == flag_t::INTERFACE_TO_EMPTY) {
                    distributeSingleMass(distributions, mass, massChange, length, coord,
                                         flags, fluidFraction, update_t::EMPTIED);
                    flags[curFlag] = flag_t::EMPTY;
                }
            }
        }
    }

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < mass.size(); ++i) {
        const double curDensity = computeDensity(&distributions[i * Q]);
        mass[i] += massChange[i];
        fluidFraction[i] = mass[i] / curDensity;
    }
}