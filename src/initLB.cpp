
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

void readParameters(coord_t &length, coord_t &procs, double &tau, boundary_t &boundaryConditions,
                    int &timesteps, int &timestepsPerPlotting, std::string &scenario,
                    char *parameterFile, bool verbose) {
    auto config = ConfigParser(std::string(parameterFile), verbose);

    length[0] = config.parse<int>("xlength");
    length[1] = config.parse<int>("ylength");
    length[2] = config.parse<int>("zlength");

    procs[0] = config.parse<int>("xprocs");
    procs[1] = config.parse<int>("yprocs");
    procs[2] = config.parse<int>("zprocs");

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

void initialiseFlagField(std::vector<flag_t> &flagField, const char *geometryFile,
                         const coord_t &length, const coord_t &offset, const coord_t coord,
                         const coord_t procs, boundary_t &boundaryConditions, bool verbose) {
    // TODO: Remove hardcoded dam break thingy and read from file.
    // Scenario: Dam break. A third of the size of the geometry is filled with water, the other one
    // not.
    // Note: The water isn't moving yet, might want to initialise differently!
    // TODO: Maybe initialise with velocity for dam break.
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = length[0] / 4; x < length[0] + 1; ++x) {
                flagField[indexForCell(x, y, z, length)] = flag_t::EMPTY;
            }
        }
    }

    // We surround our entire geometry with a real boundary layer.
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            setBoundaryFlag(flagField, 0, y, z, length, boundaryConditions.bcLeft);
            setBoundaryFlag(flagField, length[0] + 1, y, z, length, boundaryConditions.bcRight);
        }
    }

    for (int y = 0; y < length[1] + 2; ++y) {
        for (int x = 0; x < length[0] + 2; ++x) {
            setBoundaryFlag(flagField, x, y, 0, length, boundaryConditions.bcFront);
            setBoundaryFlag(flagField, x, y, length[2] + 1, length, boundaryConditions.bcBack);
        }
    }

    for (int z = 0; z < length[2] + 2; ++z) {
        for (int x = 0; x < length[0] + 2; ++x) {
            setBoundaryFlag(flagField, x, length[1] + 1, z, length, boundaryConditions.bcTop);
            setBoundaryFlag(flagField, x, 0, z, length, boundaryConditions.bcBottom);
        }
    }
    // Of course, most of these aren't real boundaries but should rather be
    // ghost cells.
    // That's why we replace them now with real cells, if necessary.
    for (int dim = 0; dim < 3; ++dim) {
        int curCoord[3] = {0, 0, 0};
        const int d1 = otherDims[dim][0];
        const int d2 = otherDims[dim][1];
        for (curCoord[d1] = 0; curCoord[d1] < length[d1] + 2; ++(curCoord[d1])) {
            for (curCoord[d2] = 0; curCoord[d2] < length[d2] + 2; ++(curCoord[d2])) {
                if (coord[dim] != 0) {
                    curCoord[dim] = 0;
                    setBoundaryFlag(flagField, curCoord[0], curCoord[1], curCoord[2], length,
                                    flag_t::PARALLEL_BOUNDARY);
                }
                if (coord[dim] != procs[dim] - 1) {
                    curCoord[dim] = length[dim] + 1;
                    setBoundaryFlag(flagField, curCoord[0], curCoord[1], curCoord[2], length,
                                    flag_t::PARALLEL_BOUNDARY);
                }
            }
        }
    }
}

void initialiseInterface(std::vector<double> &distributions, std::vector<double> &mass,
                         std::vector<flag_t> &flags, const coord_t &length) {
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
    flagReinit(distributions, mass, flags, filled, emptied, length);
}

std::vector<double> initialiseMassField(std::vector<flag_t> &flags, const coord_t &length) {
    auto mass = std::vector<double>(flags.size());

    // Set mass for empty cells to zero, for fluid cells to the density.
    // Interface cells are generated separately so need no special case.
    // Boundary cells are treated as empty cells, the mass shouldn't matter anyway.
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
