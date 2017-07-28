#include "freeSurface.hpp"
#include <iostream>
#include <numeric>

std::array<double, 3> computeSurfaceNormal(const std::vector<double> &distributions,
                                           const std::vector<double> &density,
                                           const coord_t &position, const coord_t &length,
                                           const std::vector<double> &mass) {
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

        const double plusDensity = density[nextNeighbour];
        const double minusDensity = density[prevNeighbour];

        double plusFluidFraction = 0.0;
        if (mass[nextNeighbour] != 0.0) {
            plusFluidFraction = plusDensity / mass[nextNeighbour];
        }

        double minusFluidFraction = 0.0;
        if (mass[prevNeighbour] != 0.0) {
            minusFluidFraction = minusDensity / mass[prevNeighbour];
        }
        normal[dim] = 0.5 * (minusFluidFraction - plusFluidFraction);
    }
    return normal;
}

void streamMass(const std::vector<double> &distributions, const std::vector<double> &density,
                const std::vector<flag_t> &flags, const coord_t &length,
                std::vector<double> &mass) {
#pragma omp parallel for
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int x = 0; x < length[0] + 2; ++x) {
                double deltaMass = 0.0;

                const coord_t curCell = coord_t{x, y, z};
                const int flagIndex = indexForCell(x, y, z, length);
                // We only consider the mass going into interface cells.
                // Empty cells have zero mass, full cells have mass 1.
                if (flags[flagIndex] != flag_t::INTERFACE)
                    continue;

                const int fieldIndex = flagIndex * Q;
                const double curDensity = density[flagIndex];
                const double curFluidFraction = mass[flagIndex] / curDensity;

                for (int i = 0; i < Q; ++i) {
                    const auto &vel = LATTICEVELOCITIES[i];
                    const auto neighCell = coord_t{x + vel[0], y + vel[1], z + vel[2]};
                    const auto neighFlag = indexForCell(neighCell, length);
                    const auto neighField = neighFlag * Q;

                    if (flags[neighFlag] == flag_t::FLUID) {
                        // Exchange interface and fluid at x + \Delta t e_i (eq. 4.2)
                        deltaMass += distributions[neighField + inverseVelocityIndex(i)] -
                                     distributions[fieldIndex + i];
                    } else if (flags[neighFlag] == flag_t::INTERFACE) {
                        const double neighDensity = density[neighFlag];
                        const double neighFluidFraction = mass[neighFlag] / neighDensity;
                        // Exchange interface and interface at x + \Delta t e_i (eq. 4.2)
                        // TODO: (maybe) substitute s_e with values from table 4.1
//                        const double s_e = distributions[neighFlag * Q + inverseVelocityIndex(i)] -
                                           distributions[fieldIndex + i];

                        const double s_e = calculateSE(distributions, flags, curCell, length, i);
                        deltaMass += s_e * 0.5 * (curFluidFraction + neighFluidFraction);
                    }
                }
                mass[flagIndex] += deltaMass; // (eq. 4.4)
            }
        }
    }
}

double calculateSE(const std::vector<double> &distributions, const std::vector<flag_t> &flags,
                   const coord_t &curCell, const coord_t &length, const int curFiIndex) {

    // TODO: dont need to check the type every time, once is enough
    //check cell type of x
    bool x_hasNoFluidNeigh = true;
    bool x_hasNoEmptyNeigh = true;
    for (int i = 0; i < Q; ++i) {
        const auto &vel = LATTICEVELOCITIES[i];
        const auto x_nb = coord_t{curCell[0] + vel[0], curCell[1] + vel[1], curCell[2] + vel[2]};
        const auto x_nb_flag = indexForCell(x_nb, length);

        if (flags[x_nb_flag] == flag_t::FLUID) {
            x_hasNoFluidNeigh = false;
        }
        if (flags[x_nb_flag] == flag_t::EMPTY) {
            x_hasNoEmptyNeigh = false;
        }
    }

    const auto &vel = LATTICEVELOCITIES[curFiIndex];
    const auto x_nb = coord_t{curCell[0] + vel[0], curCell[1] + vel[1], curCell[2] + vel[2]};

    //check for cell type of x_nb
    bool x_nb_hasNoFluidNeigh = true;
    bool x_nb_hasNoEmptyNeigh = true;
    for (int i = 0; i < Q; ++i) {
        const auto &v = LATTICEVELOCITIES[i];
        const auto xnb_nb = coord_t{x_nb[0] + v[0], x_nb[1] + v[1], x_nb[2] + v[2]};
        const auto xnb_nb_flag = indexForCell(xnb_nb, length);

        if (flags[xnb_nb_flag] == flag_t::FLUID) {
            x_nb_hasNoFluidNeigh = false;
        }
        if (flags[xnb_nb_flag] == flag_t::EMPTY) {
            x_nb_hasNoEmptyNeigh = false;
        }
    }

    bool x_isStandardCell = !(x_hasNoFluidNeigh && x_hasNoEmptyNeigh);
    bool x_nb_isStandardCell = !(x_nb_hasNoFluidNeigh && x_nb_hasNoEmptyNeigh);
    
    if ((x_isStandardCell && x_nb_isStandardCell) ||
        (x_hasNoFluidNeigh && x_nb_hasNoFluidNeigh) ||
        (x_hasNoEmptyNeigh && x_hasNoFluidNeigh)) {

        return distributions[indexForCell(x_nb, length) * Q + inverseVelocityIndex(curFiIndex)] -
                distributions[indexForCell(curCell, length) * Q + curFiIndex];
    }

    if ((x_isStandardCell && x_nb_hasNoFluidNeigh) ||
        (x_hasNoEmptyNeigh && x_nb_isStandardCell) ||
        (x_hasNoEmptyNeigh && x_nb_hasNoFluidNeigh)) {

        return distributions[indexForCell(x_nb, length) * Q + inverseVelocityIndex(curFiIndex)];
    }

    if ((x_isStandardCell && x_nb_hasNoEmptyNeigh) ||
        (x_hasNoFluidNeigh && x_nb_isStandardCell) ||
        (x_hasNoFluidNeigh && x_nb_hasNoEmptyNeigh)) {
        return -distributions[indexForCell(curCell, length) * Q + curFiIndex];
    }

    return distributions[indexForCell(x_nb, length) * Q + inverseVelocityIndex(curFiIndex)] -
           distributions[indexForCell(curCell, length) * Q + curFiIndex];
}

void getPotentialUpdates(const std::vector<double> &mass, const std::vector<double> &density,
                         gridSet_t &filled, gridSet_t &emptied, const coord_t &length) {
    // Check whether we have to convert the interface to an emptied or fluid cell.
    // Doesn't actually update the flags but pushes them to a queue.
    // We do this here so we do not have to calculate the density again.

    // Offset avoids periodically switching between filled and empty status.
    const double offset = 10e-3;
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int x = 0; x < length[0] + 2; ++x) {
                auto coord = coord_t{x, y, z};
                const int flagIndex = indexForCell(coord, length);
                // Eq. 4.7
                if (mass[flagIndex] > (1 + offset) * density[flagIndex]) {
                    filled.insert(coord);
                } else if (mass[flagIndex] < -offset * density[flagIndex]) {
                    // Emptied
                    emptied.insert(coord);
                }
            }
        }
    }
}

void interpolateEmptyCell(std::vector<double> &distributions, std::vector<double> &density,
                          const std::vector<coord_t> &toBalance, const coord_t &length,
                          const std::vector<flag_t> &flags) {
// Note: We only interpolate cells that are not emptied cells themselves!
#pragma omp parallel for
    for (size_t i = 0; i < toBalance.size(); ++i) {
        const auto &cell = toBalance[i];
        const int flagIndex = indexForCell(cell[0], cell[1], cell[2], length);
        const int cellIndex = flagIndex * Q;

        int numNeighs = 0;
        double avgDensity = 0.0;
        auto avgVel = std::array<double, 3>();
        for (int i = 0; i < Q; ++i) {
            const auto &vel = LATTICEVELOCITIES[i];
            const auto neigh = coord_t{cell[0] + vel[0], cell[1] + vel[1], cell[2] + vel[2]};

            const int neighFlagIndex = indexForCell(neigh, length);
            if (flags[neighFlagIndex] == flag_t::FLUID ||
                flags[neighFlagIndex] == flag_t::INTERFACE) {
                const int neighDistrIndex = neighFlagIndex * Q;
                const double neighDensity = density[neighFlagIndex];
                std::array<double, 3> neighVelocity;
                computeVelocity(&distributions[neighDistrIndex], neighDensity,
                                neighVelocity.data());

                ++numNeighs;
                avgDensity += neighDensity;
                avgVel[0] += neighVelocity[0];
                avgVel[1] += neighVelocity[1];
                avgVel[2] += neighVelocity[2];
            }
        }

        // Every empty cell has at least one interface cell as neighbour, otherwise we have a
        // worse problem than division by zero.
        assert(numNeighs != 0);
        avgDensity /= numNeighs;
        avgVel[0] /= numNeighs;
        avgVel[1] /= numNeighs;
        avgVel[2] /= numNeighs;

        density[flagIndex] = avgDensity; // Density of new cell is changed!
        // Note: This writes the equilibrium distribution directly into the distributions array.
        computeFeq(avgDensity, avgVel.data(), &distributions[cellIndex]);
    }
}

void flagReinit(std::vector<double> distributions, std::vector<double> &mass,
                std::vector<double> &density, gridSet_t &filled, gridSet_t &emptied,
                const coord_t &length, std::vector<flag_t> &flags) {
    // First consider all filled cells.

    // Store all new fluid cells with no valid distributions.
    // We need to process them after the cells have been converted to interface cells because
    // otherwise the interpolation
    // depends on the processing order.
    // TODO: Do we need to do this? Not sure...
    auto toBalance = std::vector<coord_t>();

    for (const auto &elem : filled) {
        // Find all neighbours of this cell.
        for (const auto &vel : LATTICEVELOCITIES) {
            coord_t neighbor = elem;
            neighbor[0] += vel[0];
            neighbor[1] += vel[1];
            neighbor[2] += vel[2];
            const int neighFlag = indexForCell(neighbor[0], neighbor[1], neighbor[2], length);
            // This neighbor is converted to an interface cell iff. it is an empty cell or a cell
            // that would become an emptied cell.
            // We need to remove it from the emptied set, otherwise we might have holes in the
            // interface.
            if (flags[neighFlag] == flag_t::EMPTY) {
                flags[neighFlag] = flag_t::INTERFACE;
                // Notice that the new interface cells don't have any valid distributions.
                // They are initialised with f^{eq}_i (p_{avg}, v_{avg}), which are the average
                // density and velocity of all neighbouring fluid and interface cells.
                toBalance.emplace_back(neighbor);
            } else if (flags[neighFlag] == flag_t::INTERFACE) {
                // Already is an interface but should not be converted to an empty cell later.
                emptied.erase(neighbor);
            }
        }
        // Now we can convert the cell itself to a fluid cell.
        const int curFlag = indexForCell(elem, length);
        flags[curFlag] = flag_t::FLUID;
    }

    // Secondly we consider all emptied cells that are not needed as interface cells.
    for (auto &elem : emptied) {
        // Convert all neighbours to interface cells.
        for (const auto &vel : LATTICEVELOCITIES) {
            coord_t neighbor = {elem[0] + vel[0], elem[1] + vel[1], elem[2] + vel[2]};
            const auto neighFlag = indexForCell(neighbor, length);
            if (flags[neighFlag] == flag_t::FLUID) {
                flags[neighFlag] = flag_t::INTERFACE;
                // We can reuse the distributions as they are still valid.
            }
        }
        // Again, cell should be marked as empty, finally.
        const auto curFlag = indexForCell(elem, length);
        flags[curFlag] = flag_t::EMPTY;
    }

    // Now we can interpolate the distributions for the new interface cells that have been former
    // empty cells.
    interpolateEmptyCell(distributions, density, toBalance, length, flags);
}

// TODO: Find better name for this!
enum class update_t { FILLED, EMPTIED };

void distributeSingleMass(const std::vector<double> &distributions, std::vector<double> &mass,
                          const std::vector<double> &density, const update_t &type,
                          const coord_t &length, const coord_t &coord, std::vector<flag_t> &flags) {
    // First determine how much mass needs to be redistributed and fix mass of converted cell.
    const int flagIndex = indexForCell(coord, length);
    double excessMass;
    if (type == update_t::FILLED) {
        // Interface -> Full cell, filled cells have mass and should have a mass equal to their
        // density.
        excessMass = mass[flagIndex] - density[flagIndex];
        assert(excessMass >= 0.0);
        mass[flagIndex] = density[flagIndex];
    } else {
        // Interface -> Empty cell, empty cells should not have any mass so all mass is excess mass.
        excessMass = mass[flagIndex];
        assert(excessMass < 0.0); // Follows from the offset!
        mass[flagIndex] = 0.0;
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
    const auto normal = computeSurfaceNormal(distributions, density, coord, length, mass);
    std::array<double, 19> weights{};

    for (size_t i = 0; i < LATTICEVELOCITIES.size(); ++i) {
        const auto &vel = LATTICEVELOCITIES[i];
        coord_t neighbor = {coord[0] + vel[0], coord[1] + vel[1], coord[2] + vel[2]};

        const int neighFlag = indexForCell(neighbor, length);
        if (flags[neighFlag] != flag_t::INTERFACE)
            continue;

        const double dotProduct = normal[0] * vel[0] + normal[1] * vel[1] + normal[2] * vel[2];
        if (type == update_t::FILLED) {
            weights[i] = std::max(0.0, dotProduct);
        } else { // EMPTIED
            weights[i] = -std::min(0.0, dotProduct);
        }
        assert(weights[i] >= 0.0);
    }

    // Step 2: Calculate normalizer (otherwise sum of weights != 1.0)
    const double normalizer =
        std::accumulate(weights.begin(), weights.end(), 0.0, std::plus<double>());

    if (normalizer == 0.0)
        return; // TODO: Is this a good thing to do? This corner case isn't mentioned in the paper.

    // Step 3: Redistribute weights. As non-interface cells have weight 0, we can just loop through
    // all cells.
    for (size_t i = 0; i < LATTICEVELOCITIES.size(); ++i) {
        const auto &vel = LATTICEVELOCITIES[i];
        coord_t neighbor = {coord[0] + vel[0], coord[1] + vel[1], coord[2] + vel[2]};

        const int neighFlag = indexForCell(neighbor, length);
        mass[neighFlag] += (weights[i] / normalizer) * excessMass;
    }
}

void distributeMass(const std::vector<double> &distributions, std::vector<double> &mass,
                    const std::vector<double> &density, gridSet_t &filled, gridSet_t &emptied,
                    const coord_t &length, std::vector<flag_t> &flags) {
    // Here we redistribute the excess mass of the cells.
    // It is important that we get a copy of the filled/emptied where all converted cells are stored
    // and no other cells.
    // This excludes emptied cells that are used as interface cells instead!

    for (auto &&coord : filled) {
        distributeSingleMass(distributions, mass, density, update_t::FILLED, length, coord, flags);
    }
    for (auto &&coord : emptied) {
        distributeSingleMass(distributions, mass, density, update_t::EMPTIED, length, coord, flags);
    }
}