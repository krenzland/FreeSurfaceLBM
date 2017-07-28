#include "VtkWriter.hpp"
#include "boundary.hpp"
#include "collision.hpp"
#include "freeSurface.hpp"
#include "initLB.hpp"
#include "streaming.hpp"
#include "timeStep.hpp"
#include <cassert>
#include <chrono>

int main(int argc, char *argv[]) {
    coord_t length;
    int timesteps;
    int timestepsPerPlotting;
    double stepSize = 1.0;
    double tau;
    std::array<double, 3> gravity;
    boundary_t boundaryConditions;
    std::string scenario;

    const auto startTime = std::chrono::high_resolution_clock::now();

    // Only show output for first rank.
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
    auto flagField = std::vector<flag_t>(num_cells);
    auto density = std::vector<double>(num_cells, 1.0);

    initialiseCollideAndStreamFields(collideField, streamField);
    initialiseFlagField(flagField, scenario, length, boundaryConditions, verbose);
    auto mass = initialiseMassField(flagField, length);
    initialiseInterface(streamField, mass, density, length, flagField);

    auto writer = VtkWriter("results/output", length);
    writer.write(collideField, mass, density, flagField, 0);

    // We only adapt the time step every few iterations. This here is a heuristic that seems to work
    // in practice.
    const int rescaleDelay =
        static_cast<int>(4.0 * std::pow(length[0] * length[1] * length[2], 1.0 / 3.0))*0 + 1;
    for (int t = 1; t < timesteps; ++t) {
        auto filled = gridSet_t();
        auto emptied = gridSet_t();

        doStreaming(collideField, streamField, mass, density, length, flagField);
        streamMass(streamField, density, flagField, length, mass); // Maybe do after normal streaming?
        std::swap(collideField, streamField);
        doCollision(collideField, mass, density, flagField, tau, gravity, length, filled, emptied);
        getPotentialUpdates(mass, density, filled, emptied, length);
        flagReinit(collideField, mass, density, filled, emptied, length, flagField);
        distributeMass(collideField, mass, density, filled, emptied, length, flagField);

        if (t % rescaleDelay == 0) {
            std::tie(tau, stepSize) =
                adaptTimestep(collideField, density, mass, flagField, tau, stepSize, gravity);
            std::cout << "It = " << t << " stepSize = " << stepSize << " tau = " << tau
                      << std::endl;
        }

        treatBoundary(collideField, flagField, boundaryConditions, length);

        if (!(t % timestepsPerPlotting)) {
            writer.write(collideField, mass, density, flagField, t);
        }
    }

    const auto endTime = std::chrono::high_resolution_clock::now();
    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    const double MLUPS = (length[0] * length[1] * length[2] * timesteps) / (elapsedTime * 10e6);
    std::cout << "Time elapsed " << elapsedTime << '\n' << "MLUPS " << MLUPS << std::endl;

    return 0;
}
