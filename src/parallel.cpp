#include "parallel.hpp"
#include "LBMHelper.hpp"
#include <assert.h>
#include <iostream>

void initializeMPI(int &rank, int &noRanks, int &argc, char ***argv) {
    MPI_Init(&argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &noRanks);
}

MPI_Comm initializeGrid(int rank, int noRanks, const coord_t &procs, coord_t &coords,
                        bool verbose) {
    MPI_Comm mpiGrid;

    const int periodicDim[3] = {false, false, false}; // No periodic grid!
    const bool reorderTopology =
        false; // Do not reorder the ranks, this might lead to confusing stuff.
    MPI_Barrier(MPI_COMM_WORLD);
    if (procs[0] * procs[1] * procs[2] != noRanks) {
        if (verbose) {
            printf("Creation of Cartesian grid failed. You tried to use (%d x "
                   "%d x %d) processes, "
                   "but there are %d MPI-ranks.\n",
                   procs[0], procs[1], procs[2], noRanks);
        }
        exit(-1);
    }
    MPI_Cart_create(MPI_COMM_WORLD, 3, procs.data(), periodicDim, reorderTopology, &mpiGrid);
    MPI_Cart_coords(mpiGrid, rank, 3, coords.data());
    return mpiGrid;
}
std::pair<parallelBuffers_t, parallelBuffers_t>
initializeCommunicationBuffers(const coord_t &proc, const coord_t &length,
                               const std::array<int, 6> &neighborsIdx) {
    parallelBuffers_t bufferIn, bufferOut;
    for (int i = 0; i < 6; ++i) {
        // We only need to allocate a buffer when there actually is a neighbour!
        if (neighborsIdx[i] != MPI_PROC_NULL) {
            const auto &dims = otherDims[i / 2];
            const int bufferLen =
                5 * (length[dims[0]] / proc[dims[0]] + 2) * (length[dims[1]] / proc[dims[1]] + 2);

            bufferIn[i] = std::vector<double>(bufferLen);
            bufferOut[i] = std::vector<double>(bufferLen);
        } else {
            // Otherwise use sane defaults.
            bufferIn[i] = std::vector<double>();
            bufferOut[i] = std::vector<double>();
        }
    }
    return std::pair<parallelBuffers_t, parallelBuffers_t>(std::move(bufferIn),
                                                           std::move(bufferOut));
}

void decomposeDomain(coord_t &length, const coord_t &procs, bool verbose) {
    for (size_t i = 0; i < length.size(); ++i) {
        if (length[i] % procs[i] == 0) {
            // Processes divide domain in equal parts!
            length[i] /= procs[i];
        } else {
            if (verbose) {
                std::cout << "Domain is not divided equally by processes. Make "
                             "sure that length is an integer multiple of "
                             "the number of processes.\n!";
            }
            exit(-1);
        }
    }
}

int getNeighborMPI(MPI_Comm comm, const coord_t coords, const coord_t &procs, int dim, int offset) {
    coord_t newCoords;

    for (int i = 0; i < 3; ++i) {
        int newCoord = coords[i];
        if (i == dim)
            newCoord += offset;
        if (newCoord < 0 || newCoord >= procs[i]) {
            // For all non valid coordinates return a value indicating that we
            // have no neighbour.
            return MPI_PROC_NULL;
        }
        newCoords[i] = newCoord;
    }
    int neighborRank;
    MPI_Cart_rank(comm, newCoords.data(), &neighborRank);

    return neighborRank;
}

// Each entry: sign (e.g. left or right?) , dimension \in {0,1,2}
const static int boundaryPos[6][2] = {{-1, 0}, {+1, 0}, {-1, 1}, {+1, 1}, {-1, 2}, {+1, 2}};

std::array<int, 6> getNeighborsMPI(MPI_Comm comm, const coord_t &coords, const coord_t &procs) {
    std::array<int, 6> neighbors;
    for (size_t i = 0; i < neighbors.size(); ++i) {
        neighbors[i] = getNeighborMPI(comm, coords, procs, boundaryPos[i][1], boundaryPos[i][0]);
    }
    return neighbors;
}
const static int relevantDistributions[6][5] = {
    {1, 5, 8, 11, 15},    // right
    {3, 7, 10, 13, 17},   // left
    {0, 5, 6, 7, 14},     // down
    {4, 11, 12, 13, 18},  // up
    {0, 1, 2, 3, 4},      // back
    {14, 15, 16, 17, 18}, // forth
};

const static int reverseDistribution[6] = {1, 0, 3, 2, 5, 4};

enum class direction_t { SEND, RECEIVE };

void synchBuffer(int idx, std::vector<double> &collideField, const std::vector<flag_t> &flagField,
                 std::vector<double> &buffer, const coord_t &length, direction_t direction) {
    const int noDistr = 5;
    const int dim = boundaryPos[idx][1];
    assert(dim >= 0 && dim <= 2);
    const int firstDim = otherDims[dim][0];
    const int secondDim = otherDims[dim][1];

    int curCoord[3] = {0, 0, 0};
    // If we inject stuff the boundary is on the opposite site!
    // We do not need to manually reverse it, this is taken care of by the
    // boundaryPos array.
    bool isFirst = boundaryPos[idx][0] == -1;

    if (direction == direction_t::SEND) {
        // We stream from inner cells directly, not from the parallel
        // boundaries.
        curCoord[dim] = isFirst ? 1 : (length[dim]);
    } else {
        // But we insert into parallel boundary cells.
        curCoord[dim] = isFirst ? 0 : (length[dim] + 1);
    }

    int bufferPos = 0; // Tracks where we are in the buffer
    // Iterate over the entire boundary wall. Order is important, otherwise we
    // might flip things around.
    for (curCoord[firstDim] = 0; curCoord[firstDim] < length[firstDim] + 2; ++curCoord[firstDim]) {
        for (curCoord[secondDim] = 0; curCoord[secondDim] < length[secondDim] + 2;
             ++curCoord[secondDim]) {
            const int flagIdx = indexForCell(curCoord[0], curCoord[1], curCoord[2], length);
            for (int k = 0; k < noDistr; ++k) {
                // The distributions in the buffer are continuously ordered.
                // That is not the case for the collide field.
                if (direction == direction_t::SEND) {
                    const int curDistr = relevantDistributions[idx][k];

                    buffer[bufferPos + k] = collideField[flagIdx * Q + curDistr];
                } else {
                    assert(flagField[flagIdx] != flag_t::FLUID);
                    const int curDistr = relevantDistributions[reverseDistribution[idx]][k];

                    collideField[flagIdx * Q + curDistr] = buffer[bufferPos + k];
                }
            }
            bufferPos += noDistr;
        }
    }
}

void extract(int idx, std::vector<double> &collideField, const std::vector<flag_t> &flagField,
             std::vector<double> &bufferOut, coord_t &length) {
    synchBuffer(idx, collideField, flagField, bufferOut, length, direction_t::SEND);
}
void inject(int idx, std::vector<double> &collideField, const std::vector<flag_t> &flagField,
            std::vector<double> &bufferIn, const coord_t &length) {
    synchBuffer(idx, collideField, flagField, bufferIn, length, direction_t::RECEIVE);
}
void swap(std::vector<double> &bufferIn, std::vector<double> &bufferOut,
          const std::array<int, 6> &neighbors, int neighborIdx, MPI_Comm comm) {
    MPI_Sendrecv(bufferOut.data(), bufferOut.size(), MPI_DOUBLE, neighbors[neighborIdx], 0,
                 bufferIn.data(), bufferIn.size(), MPI_DOUBLE, neighbors[neighborIdx], MPI_ANY_TAG,
                 comm, MPI_STATUS_IGNORE);
}

void swapNonBlocking(double *bufferIn, double *bufferOut, int bufferLen, int *neighbors,
                     int neighborIdx, MPI_Comm comm, MPI_Request *req) {
    MPI_Request sendReq; // We ignore that one.
    MPI_Isend(bufferOut, bufferLen, MPI_DOUBLE, neighbors[neighborIdx], 0, comm, &sendReq);
    MPI_Irecv(bufferIn, bufferLen, MPI_DOUBLE, neighbors[neighborIdx], 0, comm, req);
}

std::pair<coord_t, coord_t> lengthWithoutGhostCells(coord_t length, const coord_t &coords,
                                                    const coord_t &procs) {
    coord_t offset;
    for (int i = 0; i < 3; ++i) {
        offset[i] = (length[i] - 1) * coords[i];

        // Both cases can be true at the same time, if there is only one process
        // in this direction!
        if (coords[i] == 0) {
            ++(length[i]);
        } else {
            ++(offset[i]); // There is a boundary left of the thread!
        }
        if (coords[i] == procs[i] - 1) {
            ++(length[i]);
        }
    }
    return std::pair<coord_t, coord_t>(std::move(length), std::move(offset));
}
