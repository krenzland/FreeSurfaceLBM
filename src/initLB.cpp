
#include "initLB.hpp"
#include "ConfigParser.hpp"
#include "LBMHelper.hpp"
#include "freeSurface.hpp"
#include "scenarios/Scenario.hpp"
#include <scenarios/CornerDamBreak.hpp>
#include <scenarios/DamBreak.hpp>
#include <scenarios/FallingDrop.hpp>
#include <scenarios/MiddleWall.hpp>
#include <scenarios/OnlyWater.hpp>
#include <scenarios/MultipleDrops.hpp>
#include <scenarios/HoleInContainer.hpp>

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

void readParameters(coord_t &length, double &tau, double &smagorinskyConstant,
                    boundary_t &boundaryConditions, int &timesteps, int &timestepsPerPlotting,
                    std::unique_ptr<Scenario> &scenario, char *parameterFile, bool verbose,
                    std::array<double, 3> &gravity) {
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
    smagorinskyConstant = config.parse<double>("smagorinskyConstant");

    auto scenarioName = config.parse<std::string>("scenario");
    if (scenarioName == "onlyWater") {
        scenario = std::make_unique<OnlyWater>();
    } else if (scenarioName == "damBreak") {
        scenario = std::make_unique<DamBreak>(config);
    } else if (scenarioName == "fallingDrop") {
        scenario = std::make_unique<FallingDrop>(config);
    } else if (scenarioName == "cornerDamBreak") {
        scenario = std::make_unique<CornerDamBreak>(config);
    } else if (scenarioName == "middleWall") {
        scenario = std::make_unique<MiddleWall>(config);
    } else if (scenarioName == "multipleDrops") {
        scenario = std::make_unique<MultipleDrops>(config);
    } else if (scenarioName == "holeInContainer") {
        scenario = std::make_unique<HoleInContainer>(config);
    }else {
        throw std::invalid_argument("Invalid scenario!");
    }

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
    if (flagField[index] != flag_t::NO_SLIP) {
        flagField[index] = value;
    }
}

void initialiseFlagField(std::vector<flag_t> &flagField, std::unique_ptr<Scenario> scenario,
                         boundary_t &boundaryConditions, bool verbose, const coord_t &length) {
    scenario->getFlagField(flagField, length);
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

void initialiseInterface(std::vector<double> &distributions, std::vector<double> &mass, std::vector<double> &fluidFraction,
                         const coord_t &length, std::vector<flag_t> &flags) {
    auto newFlags = flags;
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = 1; x < length[0] + 1; ++x) {
                const auto coord = coord_t{x, y, z};
                const int flagIndex = indexForCell(coord, length);
                if (flags[flagIndex] == flag_t::FLUID) {
                    for (auto &vel : LATTICEVELOCITIES) {
                        const coord_t neigh = coord_t{x + vel[0], y + vel[1], z + vel[2]};
                        const int neighFlag = indexForCell(neigh, length);
                        if (flags[neighFlag] == flag_t::EMPTY) {
                            newFlags[neighFlag] = flag_t::INTERFACE;
                            mass[neighFlag] = 0.5;
                            fluidFraction[neighFlag] = 0.5;
                        }
                    }
                }
            }
        }
    }
    flags = std::move(newFlags);
}

std::pair<std::vector<double>, std::vector<double>> initialiseMassAndFluidFractionFields(std::vector<flag_t> &flags,
                                                                                         const coord_t &length) {
    auto mass = std::vector<double>(flags.size());
    auto fluidFraction = std::vector<double>(flags.size());

// Set mass for empty cells to zero, for fluid cells to the density.
// Interface cells are generated separately so need no special case.
// Boundary cells are treated as empty cells, the mass shouldn't matter anyway.
#pragma omp parallel for
    for (size_t i = 0; i < flags.size(); ++i) {
        if (flags[i] == flag_t::FLUID) {
            // Density in first timestep is 1 for fluid cells.
            mass[i] = 1.0;
            fluidFraction[i] = 1.0;
        } else if (flags[i] == flag_t::INTERFACE) {
            mass[i] = 0.5; // Arbitrary value, doesn't get converted too soon.
            fluidFraction[i] = 0.5;
        } else {
            mass[i] = 0.0;
            fluidFraction[i] = 0.0;
        }
    }
    return std::pair<std::vector<double>, std::vector<double>>(std::move(mass), std::move(fluidFraction));
}
