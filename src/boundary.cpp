#include "boundary.hpp"
#include "LBMHelper.hpp"
#include "computeCellValues.hpp"
#include <cassert>

struct mirror_t {
    int from; // Index of fluid distribution
    int to;   // Index of boundary distribution
};

/*
 * These are all velocities that have to be mirrored so that the tangential
 * velocity is preserved
 * and s.t. the normal velocity vanishes.
 * The first index corresponds to the position of the fluid cell w.r.t the
 * boundary, the second index
 * maps the distributions.
 * Note that for +/- positions this map is symmetric.
 * */
static const mirror_t mirroredIdx[6][5] = {
    {{3, 1}, {7, 5}, {10, 8}, {13, 11}, {17, 15}}, // x+
    {{4, 0}, {11, 5}, {12, 6}, {13, 7}, {18, 14}}, // y+
    {{14, 0}, {15, 1}, {16, 2}, {17, 3}, {18, 4}}, // z+
    {{1, 3}, {5, 7}, {8, 10}, {11, 13}, {15, 17}}, // x-
    {{0, 4}, {5, 11}, {6, 12}, {7, 13}, {14, 18}}, // y-
    {{0, 14}, {1, 15}, {2, 16}, {3, 17}, {4, 18}}, // z-
};

void mirrorFreeSlip(std::vector<double> &collideField, const std::vector<flag_t> &flagField, int x,
                    int y, int z, int offset_x, int offset_y, int offset_z, const coord_t &length) {
    const int flagIndex = indexForCell(x, y, z, length);
    const int collideIndex = flagIndex * Q;

    // We only consider non-diagonal direct neighbours here.
    assert((abs(offset_x) + abs(offset_y) + abs(offset_z)) == 1);

    // Check whether this neighbour exists and is a fluid.
    // Otherwise we don't have to consider this cell.
    const int nX = x + offset_x;
    const int nY = y + offset_y;
    const int nZ = z + offset_z;
    if (nX < 0 || nY < 0 || nZ < 0 || nX > length[0] + 1 || nY > length[1] + 1 ||
        nZ > length[2] + 1)
        return;

    const int neighborFlag = indexForCell(nX, nY, nZ, length);
    if (flagField[neighborFlag] != flag_t::FLUID)
        return;

    // We now know that we have a relevant fluid cell as a neighbour.
    const int neighborCollide = neighborFlag * Q;

    // This holds because we know that all offsets add up to either 1 or -1.
    int mirrorOffset = abs(offset_y) + abs(offset_z) * 2;
    if ((offset_x == -1) || (offset_y == -1) || (offset_z == -1)) {
        mirrorOffset += 3;
    }
    assert(mirrorOffset >= 0 && mirrorOffset <= 5);

    // Now we only have to switch distributions according to our table.
    double *boundaryCell = &collideField[collideIndex];
    double *fluidCell = &collideField[neighborCollide];
    const mirror_t *mirror = mirroredIdx[mirrorOffset];

    for (int i = 0; i < 5; ++i) {
        boundaryCell[mirror[i].from] = fluidCell[mirror[i].to];
    }
}

void singleFreeSlip(std::vector<double> &collideField, const std::vector<flag_t> &flagField, int x,
                    int y, int z, const coord_t &length) {
    // We have to consider six possible fluid cells, with offsets of -1/1 in
    // {x,y,z}-directions.
    mirrorFreeSlip(collideField, flagField, x, y, z, 1, 0, 0, length);
    mirrorFreeSlip(collideField, flagField, x, y, z, -1, 0, 0, length);
    mirrorFreeSlip(collideField, flagField, x, y, z, 0, 1, 0, length);
    mirrorFreeSlip(collideField, flagField, x, y, z, 0, -1, 0, length);
    mirrorFreeSlip(collideField, flagField, x, y, z, 0, 0, 1, length);
    mirrorFreeSlip(collideField, flagField, x, y, z, 0, 0, -1, length);
}

void singleNoSlip(std::vector<double> &collideField, const std::vector<flag_t> &flagField,
                  double const *wallVelocity, int x, int y, int z, const coord_t &length) {
    const int flagIndex = indexForCell(x, y, z, length);
    assert(flagField[flagIndex] != flag_t::FLUID);
    const int collideIndex = flagIndex * Q;
    // We loop over all possible velocities.
    // If they point into a cell that is not accessible, ignore this velocity.
    // Otherwise, we update it to reflect the boundary condition.
    for (int i = 0; i < Q; ++i) {
        const auto &vel = LATTICEVELOCITIES[i];
        const int nX = x + vel[0];
        const int nY = y + vel[1];
        const int nZ = z + vel[2];
        if (nX <= 0 || nY <= 0 || nZ <= 0 || nX > length[0] || nY > length[1] || nZ > length[2]) {
            // There is no neighbour in that direction, velocity doesn't matter
            // then.
            continue;
        }

        // Calculate the position of the cell the velocity points to.
        const int neighFlagIndex = indexForCell(nX, nY, nZ, length);
        const int neighCollideIndex = neighFlagIndex * Q;

        // This is a relevant direction.
        // First check, whether it's a wall, or not.
        double acceleration = 0.0;
        if (flagField[flagIndex] == flag_t::MOVING_WALL) {
            const double density = computeDensity(&collideField[neighCollideIndex]);
            const double adjusted_velocity =
                (wallVelocity[0] * vel[0] + wallVelocity[1] * vel[1] + wallVelocity[2] * vel[2]) /
                (C_S * C_S);
            acceleration = 2 * LATTICEWEIGHTS[i] * density * adjusted_velocity;
        }
        // We now reflect the velocity.
        // To do that, we first find the discrete vel. that points in the
        // opposite direction.
        const int inv = inverseVelocityIndex(i);
        assert(collideIndex % Q == 0 &&
               neighCollideIndex % Q == 0); // Make sure they point to the start of a fluid cell.
        collideField[collideIndex + i] = collideField[neighCollideIndex + inv] + acceleration;
        if (flagField[neighFlagIndex] == flag_t::EMPTY) {
            collideField[collideIndex + i] = 0.0;
        }
    }
}

void singleInflow(double &boundaryCell, double const *velocityIn) {
    const double pressureRef = 1.0;
    computeFeq(pressureRef, velocityIn, &boundaryCell);
}

void singleOutflow(std::vector<double> &collideField, double refDensity, int x, int y, int z,
                   const coord_t &length) {
    const int collideIndex = indexForCell(x, y, z, length) * Q;
    double feq[19] = {0.0};
    double incVel[3] = {0.0};
    double incDensity = 0.0;
    for (int i = 0; i < Q; ++i) {
        const auto &vel = LATTICEVELOCITIES[i];
        const int nX = x + vel[0];
        const int nY = y + vel[1];
        const int nZ = z + vel[2];
        if (nX <= 0 || nY <= 0 || nZ <= 0 || nX > length[0] || nY > length[1] || nZ > length[2]) {
            // There is no neighbour in that direction, velocity doesn't matter
            // then.
            continue;
        }
        const int neighCollideIndex = indexForCell(nX, nY, nZ, length) * Q;
        incDensity = computeDensity(&collideField[neighCollideIndex]);
        computeVelocity(&collideField[neighCollideIndex], incDensity, incVel);
        // We only need two values of the eq.distr., but this is easier.
        computeFeq(refDensity, incVel, feq);

        const int inv = inverseVelocityIndex(i);
        collideField[collideIndex + i] = feq[inv] + feq[i] - collideField[neighCollideIndex + inv];
    }
}

void treatBoundary(std::vector<double> &collideField, const std::vector<flag_t> &flagField,
                   const boundary_t &boundaryConditions, const coord_t &length) {
#pragma omp parallel for schedule(static)
    for (int x = 0; x < length[0] + 2; ++x) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int z = 0; z < length[2] + 2; ++z) {
                const int flagIndex = indexForCell(x, y, z, length);
                const int collideIndex = flagIndex * Q;
                const flag_t type = flagField[flagIndex];
                switch (type) {
                case flag_t::FLUID:
                case flag_t::EMPTY:
                case flag_t::INTERFACE:
                case flag_t::INTERFACE_TO_EMPTY:
                case flag_t::INTERFACE_TO_FLUID:
                    break;
                case flag_t::NO_SLIP:
                    singleNoSlip(collideField, flagField, boundaryConditions.velocityWall, x, y, z,
                                 length);
                    break;
                case flag_t::MOVING_WALL:
                    singleNoSlip(collideField, flagField, boundaryConditions.velocityWall, x, y, z,
                                 length);
                    break;
                case flag_t::INFLOW:
                    singleInflow(collideField[collideIndex], boundaryConditions.velocityIn);
                    break;
                case flag_t::OUTFLOW:
                    singleOutflow(collideField, 1.0, x, y, z, length);
                    break;
                case flag_t::PRESSURE_IN:
                    singleOutflow(collideField, boundaryConditions.pressureIn, x, y, z, length);
                    break;
                case flag_t::FREE_SLIP:
                    singleFreeSlip(collideField, flagField, x, y, z, length);
                    break;
                }
            }
        }
    }
}
