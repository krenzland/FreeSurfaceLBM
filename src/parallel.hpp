#ifndef _PARALLEL_H_
#define _PARALLEL_H_

#include "LBDefinitions.hpp"
#include "mpi.h"
#include <stdbool.h>
#include <vector>

using parallelBuffers_t = std::array<std::vector<double>, 6>;

void initializeMPI(int &rank, int &noRanks, int &argc, char ***argv);

MPI_Comm initializeGrid(int rank, int noRanks, const coord_t &procs, coord_t &coords, bool verbose);

std::pair<parallelBuffers_t, parallelBuffers_t>
initializeCommunicationBuffers(const coord_t &proc, const coord_t &length,
                               const std::array<int, 6> &neighborsIdx);

void decomposeDomain(coord_t &length, const coord_t &procs, bool verbose);

std::array<int, 6> getNeighborsMPI(MPI_Comm comm, const coord_t &coords, const coord_t &procs);

void extract(int idx, std::vector<double> &collideField, const std::vector<flag_t> &flagField,
             std::vector<double> &bufferOut, coord_t &length);

void inject(int idx, std::vector<double> &collideField, const std::vector<flag_t> &flagField,
            std::vector<double> &bufferIn, const coord_t &length);

void swap(std::vector<double> &bufferIn, std::vector<double> &bufferOut,
          const std::array<int, 6> &neighbors, int neighborIdx, MPI_Comm comm);

void swapNonBlocking(double *bufferIn, double *bufferOut, int bufferLen, int *neighbors,
                     int neighborIdx, MPI_Comm comm, MPI_Request *req);

std::pair<coord_t, coord_t> lengthWithoutGhostCells(coord_t length, const coord_t &coords,
                                                    const coord_t &procs);

#endif