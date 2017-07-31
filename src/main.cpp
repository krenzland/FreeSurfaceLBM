#include "VtkWriter.hpp"
#include "boundary.hpp"
#include "collision.hpp"
#include "freeSurface.hpp"
#include "initLB.hpp"
#include "streaming.hpp"
#include "timeStep.hpp"
#include <bits/unique_ptr.h>
#include <cassert>
#include <chrono>
#include <iomanip>

int main(int argc, char *argv[]) {
    coord_t length;
    int timesteps;
    int timestepsPerPlotting;
    double stepSize = 1.0;
    double smagorinskyConstant;
    double tau;
    std::array<double, 3> gravity;
    boundary_t boundaryConditions;
    auto scenario = std::unique_ptr<Scenario>(nullptr);

    const auto startTime = std::chrono::high_resolution_clock::now();

    const bool verbose = true;

    // Parse the scenario configuration.
    assert(argc == 2); // 1: parameter file

    readParameters(length, tau, smagorinskyConstant, boundaryConditions, timesteps,
                   timestepsPerPlotting, scenario, argv[1], verbose, gravity);

    // Surround with ghost/boundary cells.
    const size_t num_cells = (size_t)(length[0] + 2) * (length[1] + 2) * (length[2] + 2);

    // Initialise all fields.
    auto collideField = std::vector<double>(num_cells * Q);
    auto streamField = std::vector<double>(num_cells * Q);
    auto flagField = std::vector<flag_t>(num_cells, flag_t::FLUID);
    auto neighborhood = std::vector<neighborhood_t>(num_cells, neighborhood_t::STANDARD);

    initialiseCollideAndStreamFields(collideField, streamField);
    initialiseFlagField(flagField, std::move(scenario), boundaryConditions, verbose, length);
    std::vector<double> mass, fluidFraction;
    std::tie(mass, fluidFraction) = initialiseMassAndFluidFractionFields(flagField, length);
    initialiseInterface(streamField, mass, fluidFraction, length, flagField);

    auto writer = VtkWriter("results/output", length);
    writer.write(collideField, mass, flagField, stepSize, 0, fluidFraction);

    int realTimeSteps = 0;
    double lastOutput = 0.0;
    int fileNum = 1;
    // We allow no step size increase after a step size decrease!
    int increaseDelay = static_cast<int>(std::pow(flagField.size(), 1.0 / 3.0) * 4.0);
    int increaseNext = 1;
    for (double t = 1; t < timesteps; t += stepSize) {
        realTimeSteps++;

        doStreaming(collideField, streamField, mass, length, flagField, neighborhood,
                    fluidFraction);
        streamMass(collideField, fluidFraction, length, mass, neighborhood, flagField);
        std::swap(collideField, streamField);
        doCollision(collideField, mass, fluidFraction, smagorinskyConstant, gravity, length, tau, flagField);
        getPotentialUpdates(mass, fluidFraction, flagField, neighborhood, length);
        flagReinit(collideField, mass, fluidFraction, length, flagField);
        distributeMass(collideField, mass, length, flagField, fluidFraction);

        const double stepSizeBefore = stepSize;
        const bool allowIncrease = realTimeSteps > increaseNext;
        std::tie(tau, stepSize) = adaptTimestep(collideField, fluidFraction, mass, flagField, gravity, stepSize, tau,
                                                smagorinskyConstant, allowIncrease);
        if (stepSize != stepSizeBefore) {
            std::cout << "It = " << std::setprecision(6) << std::setw(10) << realTimeSteps
                      << " realTime " << std::setw(10) << t << std::setprecision(3)
                      << " changed to stepSize = " << std::setw(10) << stepSize
                      << " and tau = " << std::setw(10) << tau << std::endl;
        }
        if (stepSize < stepSizeBefore) {
            // Just decreased step size. Don't increase during the next few iterations.
            // This is done to avoid a heavily oscillating step size while preserving the safety
            // guarantees.
            increaseNext = realTimeSteps + increaseDelay;
        }

        treatBoundary(collideField, flagField, boundaryConditions, length);

        if ((t - lastOutput) > timestepsPerPlotting) {
            lastOutput = t;
            writer.write(collideField, mass, flagField, stepSize, fileNum++, fluidFraction);
        }
    }

    const auto endTime = std::chrono::high_resolution_clock::now();
    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    const double MLUPS =
        (1.0 * length[0] * length[1] * length[2] * realTimeSteps) / (elapsedTime * 10e6);
    std::cout << "(Wall) Time elapsed " << elapsedTime << "\nMLUPS " << MLUPS
              << "\n(Simulation) Time elapsed " << timesteps << "\nNumber of timesteps "
              << realTimeSteps << "\nAverage step size " << 1.0 * timesteps / realTimeSteps
              << std::endl;

    return 0;
}
