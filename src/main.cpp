#include "VtkWriter.hpp"
#include "boundary.hpp"
#include "collision.hpp"
#include "initLB.hpp"
#include "mpi.h"
#include "parallel.hpp"
#include "streaming.hpp"
#include "freeSurface.hpp"
#include <cassert>

/***
 * The flag field already contains the boundary cells, so no +2 is needed.
 * There is also no error check for the dimension. Its expected that the
 * xlengths are the same
 * as the dimensions in the geometry field
*/

int main(int argc, char *argv[]) {
    coord_t length;
    coord_t procs; // Number of processes for each direction.
    int timesteps;
    int timestepsPerPlotting;
    double tau;
    boundary_t boundaryConditions;
    std::string scenario;

    int rank, noRanks;
    initializeMPI(rank, noRanks, argc, &argv);
    assert(noRanks == 1); // TODO: Support MPI!
    std::cout << "Rank=" << rank << " of size " << noRanks << " initialised." << std::endl;
    const double startTime = MPI_Wtime();

    // Only show output for first rank.
    const bool verbose = rank == 0;

    // Parse the scenario configuration.
    assert(argc == 2); // 1: parameter file

    MPI_Barrier(MPI_COMM_WORLD); // Make output readable!
    readParameters(length, procs, tau, boundaryConditions, timesteps, timestepsPerPlotting,
                   scenario, argv[1], verbose);
    MPI_Barrier(MPI_COMM_WORLD);

    coord_t ourCoords; // The coordinates of our process in the new cartesian grid.
    // Create the cartesian grid, without reordering of ranks!
    const MPI_Comm mpiGrid = initializeGrid(rank, noRanks, procs, ourCoords, verbose);

    // We now have to find our neighbours.
    // Order of neighbors:
    // [0:left,1:right,2:top,3:bottom,4:front,5:back]
    std::array<int, 6> neighborsIdx = getNeighborsMPI(mpiGrid, ourCoords, procs);

    // We need 4 buffers per axis to communicate ghost cells with our neighbors.
    // Ordered the same way as neighborsIdx.
    // Their size depends on the size of the complete geometry, not only on the
    // local part.
    parallelBuffers_t bufferIn, bufferOut;
    std::tie(bufferIn, bufferOut) = initializeCommunicationBuffers(procs, length, neighborsIdx);
    decomposeDomain(length, procs, verbose);

    // Surround with ghost/boundary cells.
    const size_t num_cells = (size_t)(length[0] + 2) * (length[1] + 2) * (length[2] + 2);

    // Initialise all fields.
    auto collideField = std::vector<double>(num_cells * Q);
    auto streamField = std::vector<double>(num_cells * Q);
    auto flagField = std::vector<flag_t>(num_cells);

    // For visualisation and the flag field we need the length without ghost
    // cells, but with walls.
    coord_t realLength, offset;
    std::tie(realLength, offset) = lengthWithoutGhostCells(length, ourCoords, procs);

    initialiseCollideAndStreamFields(collideField, streamField);
    initialiseFlagField(flagField, scenario.c_str(), length, offset, ourCoords, procs,
                        boundaryConditions, verbose);
    auto mass = std::vector<double>(num_cells); // TODO: Extract to function and initialise properly!

    auto writer = VtkWriter("results/output", length, realLength, offset, ourCoords);
    writer.write(collideField, flagField, 0);

    for (int t = 1; t < timesteps; ++t) {
        auto filled = gridSet_t();
        auto emptied = gridSet_t();

        for (int i = 0; i < 6; ++i) {
            if (neighborsIdx[i] != MPI_PROC_NULL) {
                extract(i, collideField, flagField, bufferOut[i], length);
                swap(bufferIn[i], bufferOut[i], neighborsIdx, i, mpiGrid);
                inject(i, collideField, flagField, bufferIn[i], length);
            }
        }
        streamMass(collideField, flagField, mass, length); // Maybe do after normal streaming?
        doStreaming(collideField, streamField, flagField, length);
        std::swap(collideField, streamField);
        // TODO: Reconstruct boundaries for interface cells.
        doCollision(collideField, mass, flagField, tau, length, filled, emptied);
        flagReinit(mass, flagField, filled, emptied, length); // TODO: Finish implementation!
        treatBoundary(collideField, flagField, boundaryConditions, length);

        if (!(t % timestepsPerPlotting)) {
            writer.write(collideField, flagField, t);
        }
    }

    // Make sure that MPI has finished.
    const double endTime = MPI_Wtime();
    const double elapsedTime = endTime - startTime;
    const double MLUPS_single =
        (length[0] * length[1] * length[2] * timesteps) / (elapsedTime * 10e6);
    const double MLUPS_total = MLUPS_single * noRanks;
    if (rank == 0) {
        std::cout << "Time elapsed " << elapsedTime << '\n'
                  << "MLUPS (single process): " << MLUPS_single << "\nOn " << noRanks
                  << " processes " << MLUPS_total << " MLUPS." << std::endl;
    }

    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized)
        MPI_Finalize();

    return 0;
}
