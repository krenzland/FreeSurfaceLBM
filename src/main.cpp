#include "VtkWriter.hpp"
#include "boundary.hpp"
#include "collision.hpp"
#include "freeSurface.hpp"
#include "initLB.hpp"
#include "streaming.hpp"
#include <cassert>
#include <chrono>

int main(int argc, char *argv[]) {
    coord_t length;
    coord_t procs; // Number of processes for each direction.
    int timesteps;
    int timestepsPerPlotting;
    double tau;
    boundary_t boundaryConditions;
    std::string scenario;

    const auto startTime = std::chrono::high_resolution_clock::now();

    // Only show output for first rank.
    const bool verbose = true;

    // Parse the scenario configuration.
    assert(argc == 2); // 1: parameter file

    readParameters(length, procs, tau, boundaryConditions, timesteps, timestepsPerPlotting,
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
    writer.write(collideField, mass, flagField, 0);

    for (int t = 1; t < timesteps; ++t) {
        auto filled = gridSet_t();
        auto emptied = gridSet_t();

        streamMass(collideField, density, flagField, length,
                   mass); // Maybe do after normal streaming?
        doStreaming(collideField, streamField, mass, density, length, flagField);
        std::swap(collideField, streamField);
        // TODO: Reconstruct boundaries for interface cells.
        doCollision(collideField, mass, density, flagField, tau, length, filled, emptied);
        getPotentialUpdates(mass, density, filled, emptied, length);
        flagReinit(collideField, mass, density, filled, emptied, length,
                   flagField); // TODO: Finish implementation!
        distributeMass(collideField, mass, density, filled, emptied, length, flagField);
        treatBoundary(collideField, flagField, boundaryConditions, length);

        if (!(t % timestepsPerPlotting)) {
            writer.write(collideField, mass, flagField, t);
        }
    }

    const auto endTime = std::chrono::high_resolution_clock::now();
    const double elapsedTime =
        std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    const double MLUPS = (length[0] * length[1] * length[2] * timesteps) / (elapsedTime * 10e6);
    std::cout << "Time elapsed " << elapsedTime << '\n' << "MLUPS " << MLUPS << std::endl;

    return 0;
}
