
#include "initLB.hpp"
#include "ConfigParser.hpp"
#include "LBMHelper.hpp"
#include "freeSurface.hpp"

boundary_t readBoundaryConditions(ConfigParser &config) {
    boundary_t bc;
    bc.bcLeft = config.parse<flag_t>("bcLeft");
    bc.bcRight = config.parse<flag_t>("bcRight");
    bc.bcTop = config.parse<flag_t>("bcTop");
    bc.bcBottom = config.parse<flag_t>("bcBottom");
    bc.bcFront = config.parse<flag_t>("bcFront");
    bc.bcBack = config.parse<flag_t>("bcBack");

    bc.velocityWall[0] = config.parse<double>("velocityWallX");
    bc.velocityWall[1] = config.parse<double>("velocityWallY");
    bc.velocityWall[2] = config.parse<double>("velocityWallZ");

    bc.velocityIn[0] = config.parse<double>("velocityInX");
    bc.velocityIn[1] = config.parse<double>("velocityInY");
    bc.velocityIn[2] = config.parse<double>("velocityInZ");

    bc.pressureIn = config.parse<double>("pressureIn");
    return bc;
}

void readParameters(coord_t &length, double &tau, std::array<double, 3> &gravity,
                    boundary_t &boundaryConditions, int &timesteps, int &timestepsPerPlotting,
                    std::string &scenario, char *parameterFile, bool verbose) {
    auto config = ConfigParser(std::string(parameterFile), verbose);

    length[0] = config.parse<int>("xlength");
    length[1] = config.parse<int>("ylength");
    length[2] = config.parse<int>("zlength");

    gravity[0] = config.parse<double>("gravityX");
    gravity[1] = config.parse<double>("gravityY");
    gravity[2] = config.parse<double>("gravityZ");

    timesteps = config.parse<int>("timesteps");
    timestepsPerPlotting = config.parse<int>("timestepsPerPlotting");

    tau = config.parse<double>("tau");

    scenario = config.parse<std::string>("scenario");
    boundaryConditions = readBoundaryConditions(config);
}

void initialiseCollideAndStreamFields(std::vector<double> &collideField,
                                      std::vector<double> &streamField) {
    // set initial values
    for (size_t i = 0; i < collideField.size(); ++i) {
        size_t f_index = i % Q;
        collideField[i] = LATTICEWEIGHTS[f_index];
        streamField[i] = LATTICEWEIGHTS[f_index];
    }
}

void setBoundaryFlag(std::vector<flag_t> &flagField, int x, int y, int z, const coord_t &length,
                     flag_t value) {
    const int index = indexForCell(x, y, z, length);
    flagField[index] = value;
}

void initialiseFlagField(std::vector<flag_t> &flagField, const std::string &geometryFile,
                         const coord_t &length, boundary_t &boundaryConditions, bool verbose) {
    // TODO: Remove hardcoded dam break thingy and read from file.
    // Scenario: Dam break. A third of the size of the geometry is filled with water, the other one
    // not.
    // Note: The water isn't moving yet, might want to initialise differently!
    // TODO: Maybe initialise with velocity for dam break.

    // Breaking  dam
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = length[0] / 4; x < length[0] + 1; ++x) {
                flagField[indexForCell(x, y, z, length)] = flag_t::EMPTY;
            }
        }
    }
/* // Falling droplet
const int middleX = length[0] / 2;
const int middleY = length[1] / 2;
const int heightZ = 2 * (length[2] / 3);
const int radius = 10;
for (int z = 1; z < length[2] + 1; ++z) {
    for (int y = 1; y < length[1] + 1; ++y) {
        for (int x = 1; x < length[0] + 1; ++x) {
            flag_t entry = flag_t::EMPTY;
            const double distance = std::sqrt(
                std::pow(middleX - x, 2) + std::pow(middleY - y, 2) + std::pow(heightZ - z, 2));
            if (distance <= radius) {
                //entry = flag_t::FLUID;
            }
            if (z <= length[2]/5) {
                entry = flag_t::FLUID;
            }
            flagField[indexForCell(x, y, z, length)] = entry;
        }
    }
}
*/

// We surround our entire geometry with a real boundary layer.
#pragma omp parallel for
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            setBoundaryFlag(flagField, 0, y, z, length, boundaryConditions.bcLeft);
            setBoundaryFlag(flagField, length[0] + 1, y, z, length, boundaryConditions.bcRight);
        }
    }

#pragma omp parallel for
    for (int y = 0; y < length[1] + 2; ++y) {
        for (int x = 0; x < length[0] + 2; ++x) {
            setBoundaryFlag(flagField, x, y, 0, length, boundaryConditions.bcFront);
            setBoundaryFlag(flagField, x, y, length[2] + 1, length, boundaryConditions.bcBack);
        }
    }

#pragma omp parallel for
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int x = 0; x < length[0] + 2; ++x) {
            setBoundaryFlag(flagField, x, length[1] + 1, z, length, boundaryConditions.bcTop);
            setBoundaryFlag(flagField, x, 0, z, length, boundaryConditions.bcBottom);
        }
    }
}

void initialiseInterface(std::vector<double> &distributions, std::vector<double> &mass,
                         std::vector<double> &density, const coord_t &length,
                         std::vector<flag_t> &flags) {
    // Idea: Treat all empty cells as recently emptied cells and let the free surface code deal with
    // it.
    // While it works, it is pretty slow. It is fast enough because it happens only for the first
    // time step and is in
    // linear time to the number of empty cells.
    auto emptied = gridSet_t();
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int x = 0; x < length[0] + 2; ++x) {
                const auto coord = coord_t{x, y, z};
                const int flagIndex = indexForCell(coord, length);
                if (flags[flagIndex] == flag_t::EMPTY) {
                    flags[flagIndex] = flag_t::INTERFACE;
                    emptied.insert(coord);
                }
            }
        }
    }
    auto filled = gridSet_t(); // We have no recently filled cells.
    flagReinit(distributions, mass, density, filled, emptied, length, flags);
}

std::vector<double> initialiseMassField(std::vector<flag_t> &flags, const coord_t &length) {
    auto mass = std::vector<double>(flags.size());

// Set mass for empty cells to zero, for fluid cells to the density.
// Interface cells are generated separately so need no special case.
// Boundary cells are treated as empty cells, the mass shouldn't matter anyway.
#pragma omp parallel for
    for (size_t i = 0; i < flags.size(); ++i) {
        if (flags[i] == flag_t::FLUID) {
            // Density in first timestep is 1 for fluid cells.
            mass[i] = 1.0;
        } else if (flags[i] == flag_t::INTERFACE) {
            mass[i] = 0.5;
        } else {
            mass[i] = 0.0;
        }
    }
    return mass;
}
