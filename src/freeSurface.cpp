#include "freeSurface.hpp"

std::array<double, 3> computeSurfaceNormal(const std::vector<double> &distributions, const std::vector<double> &mass,
                                             const coord_t &position, const coord_t &length) {
    auto normal = std::array<double, 3>();
    // We approximate the surface normal element-wise using central-differences of the fluid fraction.
    for (size_t dim = 0; dim < normal.size(); ++dim) {
        auto curPosition = position;

        // We need the next and previous neighbour for dimension dim.
        curPosition[dim]++;
        const auto nextNeighbour = indexForCell(curPosition[0], curPosition[1], curPosition[2], length);
        curPosition[dim] -= 2;
        const auto prevNeighbour = indexForCell(curPosition[0], curPosition[1], curPosition[2], length);

        const double plusDensity = computeDensity(&distributions[nextNeighbour * Q]);
        const double minusDensity = computeDensity(&distributions[prevNeighbour * Q]);

        const double plusFluidFraction = plusDensity / mass[nextNeighbour];
        const double minusFluidFraction = minusDensity / mass[prevNeighbour];
        normal[dim] = 0.5 * (minusFluidFraction - plusFluidFraction);
    }
    return normal;
}

void streamMass(const std::vector<double> &distributions, const std::vector<flag_t> &flags, std::vector<double> mass,
                const coord_t &length) {
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
                        deltaMass += s_e - distributions[fieldIndex] * 0.5 * (curFluidFraction + neighFluidFraction);
                    }
                }
                mass[flagIndex] += deltaMass; // (eq. 4.4)
            }
        }
    }
}

void getPotentialUpdates(const coord_t &coord, double mass, double density, gridSet_t &filled, gridSet_t &emptied) {
    // Offset avoids periodically switching between filled and empty status.
    const double offset = 10e-3;
    // Eq. 4.7
    if (mass > (1 + offset) * density) {
        filled.insert(coord);
    } else if (mass < -offset * density){
        // Emptied
        emptied.insert(coord);
    }
}


void flagReinit(std::vector<double> &mass, std::vector<flag_t> &flags, gridSet_t &filled, gridSet_t &emptied,
                const coord_t &length) {
    // First consider all filled cells.
    while (filled.size() > 0) {
        const auto& it = filled.begin();
        const auto& elem = *it;
        // Find all neighbours of this cell.

        for (const auto& vel: LATTICEVELOCITIES) {
            coord_t neighbor = elem;
            neighbor[0] += vel[0];
            neighbor[1] += vel[1];
            neighbor[2] += vel[2];
            const auto neighFlag = indexForCell(neighbor[0], neighbor[1], neighbor[2], length);
            // This neighbor is converted to an interface cell iff. it is an empty cell or a cell that would become an
            // emptied cell.
            // We need to remove it from the emptied set then otherwise we might have holes in the interface.
            if (flags[neighFlag] == flag_t::EMPTY) {
                // No special care needed here.
                flags[neighFlag] = flag_t::INTERFACE;
            } else if (flags[neighFlag] == flag_t::INTERFACE) {
                // Already is an interface but should not be converted.
                emptied.erase(neighbor);
            }
            // Notice that the new interface cells don't have any valid distributions.
            // They are initialised with f^{eq}_i (p_{avg}, v_{avg}), which are the average density and velocity
            // of all neighbouring fluid and interface cells.
            // Note: Only interface cells that are not going to be emptied should be considered!
            // TODO: Init. interface distributions!
        }
        // Now we can convert the cell it self to a fluid cell.
        const auto curFlag = indexForCell(elem[0], elem[1], elem[2], length);
        flags[curFlag] = flag_t::FLUID;

        filled.erase(it);
    }
    // Secondly we consider all emptied cells that are not needed as interface cells.
    while (emptied.size() > 0) {
        const auto& it = filled.begin();
        const auto& elem = *it;
        for (const auto& vel: LATTICEVELOCITIES) {
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
    // TODO: Mass conservation!
}
