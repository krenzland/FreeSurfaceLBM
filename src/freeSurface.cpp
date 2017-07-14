#include "freeSurface.hpp"

std::array<double, 3> computeSurfaceNormal(const std::vector<double> &distributions,
                                           const std::vector<double> &mass, const coord_t &position,
                                           const coord_t &length) {
    auto normal = std::array<double, 3>();
    // We approximate the surface normal element-wise using central-differences of the fluid
    // fraction.
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
                    const int neighbour = neighbouring_fi_cell_index(x, y, z, i, length);

                    if (flags[neighbour] == flag_t::FLUID) {
                        // Exchange interface and fluid at x + \Delta t e_i (eq. 4.2)
                        deltaMass += distributions[neighbour * Q + inverseVelocityIndex(i)] -
                                     distributions[fieldIndex];
                    } else if (flags[neighbour] == flag_t::INTERFACE) {
                        const double neighDensity = computeDensity(&distributions[neighbour * Q]);
                        const double neighFluidFraction = mass[neighbour] / neighDensity;
                        // Exchange interface and interface at x + \Delta t e_i (eq. 4.2)
                        // TODO: (maybe) substitute s_e with values from table 4.1
                        const double s_e = distributions[neighbour * Q + inverseVelocityIndex(i)];
                        deltaMass += s_e -
                                     distributions[fieldIndex] * 0.5 *
                                         (curFluidFraction + neighFluidFraction);
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
            // worse problem than divison by zero.
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
            // that would become an
            // emptied cell.
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
    // TODO: Mass conservation!
}
