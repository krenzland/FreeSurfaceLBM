#include "freeSurface.hpp"
#include <numeric>

std::array<double, 3> computeSurfaceNormal(const std::vector<double> &distributions,
                                           const std::vector<double> &mass, const coord_t &position,
                                           const coord_t &length) {
    auto normal = std::array<double, 3>();
    // We approximate the surface normal element-wise using central-differences of the fluid
    // fraction gradient.
    for (size_t dim = 0; dim < normal.size(); ++dim) {
        auto curPosition = position;

        // We need the next and previous neighbour for dimension dim.
        curPosition[dim]++;
        const auto nextNeighbour =
            indexForCell(curPosition[0], curPosition[1], curPosition[2], length);
        curPosition[dim] -= 2;
        const auto prevNeighbour =
            indexForCell(curPosition[0], curPosition[1], curPosition[2], length);

        const double plusDensity = computeDensity(&distributions[nextNeighbour * Q]);
        const double minusDensity = computeDensity(&distributions[prevNeighbour * Q]);

        const double plusFluidFraction = plusDensity / mass[nextNeighbour];
        const double minusFluidFraction = minusDensity / mass[prevNeighbour];
        normal[dim] = 0.5 * (minusFluidFraction - plusFluidFraction);
    }
    return normal;
}

void streamMass(const std::vector<double> &distributions, const std::vector<flag_t> &flags,
                std::vector<double> mass, const coord_t &length) {
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int x = 0; x < length[0] + 2; ++x) {
                double deltaMass = 0.0;

                const int flagIndex = indexForCell(x, y, z, length);
                // We only consider the mass going into interface cells.
                // Empty cells have zero mass, full cells have mass 1.
                if (flags[flagIndex] != flag_t::INTERFACE)
                    continue;

                const int fieldIndex = flagIndex * Q;

                const double curDensity = computeDensity(&distributions[fieldIndex]);
                const double curFluidFraction = mass[flagIndex] / curDensity;
                for (int i = 0; i < Q; ++i) {
                    const auto &vel = LATTICEVELOCITIES[i];
                    const auto neighCell = coord_t{x + vel[0], y + vel[1], z + vel[2]};
                    const auto neighFlag = indexForCell(neighCell, length);

                    if (flags[neighFlag] == flag_t::FLUID) {
                        // Exchange interface and fluid at x + \Delta t e_i (eq. 4.2)
                        deltaMass += distributions[neighFlag * Q + inverseVelocityIndex(i)] -
                                     distributions[fieldIndex];
                    } else if (flags[neighFlag] == flag_t::INTERFACE) {
                        const double neighDensity = computeDensity(&distributions[neighFlag * Q]);
                        const double neighFluidFraction = mass[neighFlag] / neighDensity;
                        // Exchange interface and interface at x + \Delta t e_i (eq. 4.2)
                        // TODO: (maybe) substitute s_e with values from table 4.1
                        const double s_e = distributions[neighFlag * Q + inverseVelocityIndex(i)] -
                                distributions[fieldIndex];
                        deltaMass += s_e * 0.5 * (curFluidFraction + neighFluidFraction);
                    }
                }
                mass[flagIndex] += deltaMass; // (eq. 4.4)
            }
        }
    }
}

void getPotentialUpdates(const coord_t &coord, double mass, double density, gridSet_t &filled,
                         gridSet_t &emptied) {
    // Offset avoids periodically switching between filled and empty status.
    const double offset = 10e-3;
    // Eq. 4.7
    if (mass > (1 + offset) * density) {
        filled.insert(coord);
    } else if (mass < -offset * density) {
        // Emptied
        emptied.insert(coord);
    }
}

void interpolateEmptyCell(std::vector<double> &distributions, const std::vector<flag_t> &flags,
                          const std::vector<coord_t> &toBalance, const coord_t &length) {
    // Note: We only interpolate cells that are not emptied cells themselves!
    for (const auto &cell : toBalance) {
        const int cellIndex = indexForCell(cell[0], cell[1], cell[2], length);

        int numNeighs = 0;
        double avgDensity = 0.0;
        auto avgVel = std::array<double, 3>();
        for (int i = 0; i < Q; ++i) {
            const int neighFlagIndex =
                neighbouring_fi_cell_index(cell[0], cell[1], cell[2], i, length);
            if (flags[neighFlagIndex] == flag_t::FLUID ||
                flags[neighFlagIndex] == flag_t::INTERFACE) {
                const int neighDistrIndex = neighFlagIndex * Q;
                const double neighDensity = computeDensity(&distributions[neighDistrIndex]);
                std::array<double, 3> neighVelocity;
                computeVelocity(&distributions[neighDistrIndex], neighDensity,
                                neighVelocity.data());

                ++numNeighs;
                avgDensity += neighDensity;
                avgVel[0] += neighVelocity[0];
                avgVel[1] += neighVelocity[1];
                avgVel[2] += neighVelocity[2];
            }

            // Every empty cell has at least one interface cell as neighbour, otherwise we have a
            // worse problem than division by zero.
            assert(numNeighs != 0);
            avgDensity /= numNeighs;
            avgVel[0] /= numNeighs;
            avgVel[1] /= numNeighs;
            avgVel[2] /= numNeighs;

            // Note: This writes the equilibrium distribution directly into the distributions array.
            computeFeq(avgDensity, avgVel.data(), &distributions[cellIndex]);
        }
    }
}

void flagReinit(std::vector<double> distributions, std::vector<double> &mass,
                std::vector<flag_t> &flags, gridSet_t &filled, gridSet_t &emptied,
                const coord_t &length) {
    // First consider all filled cells.

    // Store all new fluid cells with no valid distributions.
    // We need to process them after the cells have been converted to interface cells because
    // otherwise the interpolation
    // depends on the processing order.
    // TODO: Do we need to do this? Not sure...
    auto toBalance = std::vector<coord_t>();

    while (filled.size() > 0) {
        const auto &it = filled.begin();
        const auto &elem = *it;
        // Find all neighbours of this cell.

        for (const auto &vel : LATTICEVELOCITIES) {
            coord_t neighbor = elem;
            neighbor[0] += vel[0];
            neighbor[1] += vel[1];
            neighbor[2] += vel[2];
            const auto neighFlag = indexForCell(neighbor[0], neighbor[1], neighbor[2], length);
            // This neighbor is converted to an interface cell iff. it is an empty cell or a cell
            // that would become an emptied cell.
            // We need to remove it from the emptied set, otherwise we might have holes in the
            // interface.
            if (flags[neighFlag] == flag_t::EMPTY) {
                // No special care needed here.
                flags[neighFlag] = flag_t::INTERFACE;
            } else if (flags[neighFlag] == flag_t::INTERFACE) {
                // Already is an interface but should not be converted.
                emptied.erase(neighbor);
            }
            // Notice that the new interface cells don't have any valid distributions.
            // They are initialised with f^{eq}_i (p_{avg}, v_{avg}), which are the average density
            // and velocity
            // of all neighbouring fluid and interface cells.
            // Note: Only interface cells that are not going to be emptied should be considered!
            toBalance.emplace_back(neighbor);
        }
        // Now we can convert the cell it self to a fluid cell.
        const auto curFlag = indexForCell(elem[0], elem[1], elem[2], length);
        flags[curFlag] = flag_t::FLUID;

        filled.erase(it);
    }

    // Secondly we consider all emptied cells that are not needed as interface cells.
    while (emptied.size() > 0) {
        const auto &it = filled.begin();
        const auto &elem = *it;
        // Convert all neighbours to interface cells.
        for (const auto &vel : LATTICEVELOCITIES) {
            coord_t neighbor = elem;
            neighbor[0] += vel[0];
            neighbor[1] += vel[1];
            neighbor[2] += vel[2];
            const auto neighFlag = indexForCell(neighbor[0], neighbor[1], neighbor[2], length);
            if (flags[neighFlag] == flag_t::FLUID) {
                flags[neighFlag] = flag_t::INTERFACE;
                // We can reuse the distributions as they are still valid.
            }
        }
        const auto curFlag = indexForCell(elem[0], elem[1], elem[2], length);
        flags[curFlag] = flag_t::EMPTY;

        emptied.erase(it);
    }

    // Now we can interpolate the distributions for the new interface cells that have been former
    // empty cells.
    interpolateEmptyCell(distributions, flags, toBalance, length);
}

// TODO: Find better name for this!
enum class update_t { FILLED, EMPTIED };

void distributeSingleMass(const std::vector<double> &distributions, std::vector<double> &mass,
                          std::vector<flag_t> &flags, const update_t &type, const coord_t &length,
                          const coord_t &coord) {
    // First determine how much mass needs to be redistributed and fix mass of converted cell.
    const int flagIndex = indexForCell(coord, length);
    double excessMass;
    if (type == update_t::FILLED) {
        const double density = computeDensity(&distributions[flagIndex * Q]);
        // Interface -> Full cell, filled cells have mass and should have a mass equal to their
        // density.
        excessMass = mass[flagIndex] - density;
        mass[flagIndex] = density;
    } else {
        // Interface -> Empty cell, empty cells should not have any mass so all mass is excess mass.
        excessMass = mass[flagIndex];
        mass[flagIndex] = 0.0;
    }

    /* The distribution of excess mass is surprisingly non-trivial.
       For a more detailed description, refer to pages 32f. of Thürey's thesis but here's the gist:
       We do not distribute the mass uniformly to all neighbouring interface cells but rather
       correct for balance.
       The reason for this is that the fluid interface moved beyond the current cell.
       We rebalance things by weighting the mass updates according to the direction of the interface
       normal.
       This has to be done in two steps, we first calculate all update weights, normalize them and
       then, in a second step update the weights.*/

    // Step 1: Calculate the unnormalized weights.
    const auto normal = computeSurfaceNormal(distributions, mass, coord, length);
    std::array<double, 19> weights;

    for (size_t i = 0; i < LATTICEVELOCITIES.size(); ++i) {
        const auto &vel = LATTICEVELOCITIES[i];
        coord_t neighbor = coord;
        neighbor[0] += vel[0];
        neighbor[1] += vel[1];
        neighbor[2] += vel[2];

        const int neighFlag = indexForCell(neighbor, length);
        if (flags[neighFlag] != flag_t::INTERFACE)
            continue;

        const double dotProduct = normal[0] * vel[0] + normal[1] * vel[1] + normal[2] * vel[2];
        if (type == update_t::FILLED) {
            weights[i] = dotProduct > 0 ? dotProduct : 0;
        } else {
            weights[i] = dotProduct < 0 ? -dotProduct : 0;
        }
    }

    // Step 2: Calculate normalizer (otherwise sum of weights != 1.0)
    const double normalizer =
        std::accumulate(weights.begin(), weights.end(), 0.0, std::plus<double>());

    // Step 3: Redistribute weights. As non-interface cells have weight 0, we can just loop through
    // all cells.
    for (size_t i = 0; i < LATTICEVELOCITIES.size(); ++i) {
        const auto &vel = LATTICEVELOCITIES[i];
        coord_t neighbor = coord;
        neighbor[0] += vel[0];
        neighbor[1] += vel[1];
        neighbor[2] += vel[2];

        const int neighFlag = indexForCell(neighbor, length);
        mass[neighFlag] += (weights[i] / normalizer) * excessMass;
    }
}

void distributeMass(const std::vector<double> &distributions, std::vector<double> &mass,
                    std::vector<flag_t> &flags, gridSet_t &filled, gridSet_t &emptied,
                    const coord_t &length) {
    // Here we redistribute the excess mass of the cells.
    // It is important that we get a copy of the filled/emptied where all converted cells are stored
    // and no other cells.
    // Watch out for cells that are not actually emptied but rather are new interface cells!
    // TODO: Implement this! (otherwise filled/emptied are empty sets and this method would do
    // literally nothing!)

    for (auto &&coord : filled) {
        distributeSingleMass(distributions, mass, flags, update_t::FILLED, length, coord);
    }
    for (auto &&coord : emptied) {
        distributeSingleMass(distributions, mass, flags, update_t::EMPTIED, length, coord);
    }
}