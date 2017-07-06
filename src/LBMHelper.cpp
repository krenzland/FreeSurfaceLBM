//
// Created by Andreas Reiser on 24.05.17.
//

#include "LBDefinitions.hpp"
#include "assert.h"

int indexForCell(int x, int y, int z, int *length) {
    assert(x >= 0 && x < length[0] + 2);
    assert(y >= 0 && y < length[1] + 2);
    assert(z >= 0 && z < length[2] + 2);

    return (length[0] + 2) * (length[1] + 2) * z + (length[0] + 2) * y + x;
}

int indexForCell(int x, int y, int z, const coord_t &length) {
    assert(x >= 0 && x < length[0] + 2);
    assert(y >= 0 && y < length[1] + 2);
    assert(z >= 0 && z < length[2] + 2);

    return (length[0] + 2) * (length[1] + 2) * z + (length[0] + 2) * y + x;
}

int inverseVelocityIndex(int index) {
    assert(index >= 0 && index < Q);
    return (Q - 1) - index;
}
