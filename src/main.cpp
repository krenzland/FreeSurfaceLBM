#include "VtkWriter.hpp"
#include "boundary.hpp"
#include "collision.hpp"
#include "freeSurface.hpp"
#include "initLB.hpp"
#include "streaming.hpp"
#include "timeStep.hpp"
#include <cassert>
#include <chrono>
#include <bits/unique_ptr.h>

int main(int argc, char *argv[]) {
    coord_t length;
    int timesteps;
    int timestepsPerPlotting;
    double stepSize = 1.0;
    double tau;
    std::array<double, 3> gravity;
    boundary_t boundaryConditions;
    auto scenario = std::unique_ptr<Scenario>(nullptr);

    const auto startTime = std::chrono::high_resolution_clock::now();

    const bool verbose = true;

    // Parse the scenario configuration.
    assert(argc == 2); // 1: parameter file

    readParameters(length, tau, gravity, boundaryConditions, timesteps, timestepsPerPlotting,
                   scenario, argv[1], verbose);

    // Surround with ghost/boundary cells.
    const size_t num_cells = (size_t)(length[0] + 2) * (length[1] + 2) * (length[2] + 2);

    // Initialise all fields.
    auto collideField = std::vector<double>(num_cells * Q);
    auto streamField = std::vector<double>(num_cells * Q);
    auto flagField = std::vector<flag_t>(num_cells, flag_t::FLUID);
    auto density = std::vector<double>(num_cells, 1.0);

    initialiseCollideAndStreamFields(collideField, streamField);
    initialiseFlagField(flagField, std::move(scenario), boundaryConditions,
                        verbose, length);
    auto mass = initialiseMassField(flagField, length);
    initialiseInterface(streamField, mass, density, length, flagField);
    mass = initialiseMassField(flagField, length);

    auto writer = VtkWriter("results/output", length);
    writer.write(collideField, mass, density, flagField, stepSize, stepSize);

    int realTimeSteps = 0;
    double lastOutput = 0.0;
    int fileNum = 1;
    for (double t = 1; t < timesteps; t+=stepSize) {
        realTimeSteps++;

        doStreaming(collideField, streamField, mass, density, length, flagField);
        streamMass(collideField, density, flagField, length, mass); // Maybe do after normal streaming?
        std::swap(collideField, streamField);
        doCollision(collideField, mass, density, flagField, tau, gravity, length);
        getPotentialUpdates(mass, density, length, flagField);
        flagReinit(collideField, mass, density, length, flagField);
        distributeMass(collideField, mass, density, length, flagField);

        const double stepSizeBefore = stepSize;
        // TODO: Use a delay for increasing (maybe).
        std::tie(tau, stepSize) =
                adaptTimestep(collideField, density, mass, flagField, tau, stepSize, gravity);
        if (stepSize != stepSizeBefore) {
            std::cout << "It = " << realTimeSteps << " stepSize = " << stepSize << " tau = " << tau
                      << std::endl;
        }

        treatBoundary(collideField, flagField, boundaryConditions, length);

        if ( (t-lastOutput) > timestepsPerPlotting ) {
            lastOutput = t;
            writer.write(collideField, mass, density, flagField, stepSize, fileNum++);
        }

    }

    const auto endTime = std::chrono::high_resolution_clock::now();
    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    const double MLUPS = (length[0] * length[1] * length[2] * realTimeSteps) / (elapsedTime * 10e6);
    std::cout << "(Wall) Time elapsed " << elapsedTime << "\nMLUPS " << MLUPS
              <<  "\n(Simulation) Time elapsed " << timesteps << "\nNumber of timesteps " << realTimeSteps
              << "\nAverage step size " << 1.0*timesteps/realTimeSteps << std::endl;

    return 0;
}
