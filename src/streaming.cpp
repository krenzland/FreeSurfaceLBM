#include "streaming.hpp"
#include "LBMHelper.hpp"
#include <assert.h>

int neighbouring_fi_cell_index(int x, int y, int z, int fi, const coord_t &length) {
    int new_x = x - LATTICEVELOCITIES[fi][0];
    int new_y = y - LATTICEVELOCITIES[fi][1];
    int new_z = z - LATTICEVELOCITIES[fi][2];

    return indexForCell(new_x, new_y, new_z, length);
}

void doStreaming(const std::vector<double> &collideField, std::vector<double> &streamField,
                 const std::vector<flag_t> &flagField, const coord_t &length) {
    for (int z = 0; z < length[2] + 2; ++z) {
        for (int y = 0; y < length[1] + 2; ++y) {
            for (int x = 0; x < length[0] + 2; ++x) {
                const int flag_index = indexForCell(x, y, z, length);
                if (flagField[flag_index] != flag_t::FLUID)
                    continue;

                const int field_index = flag_index * Q;
                for (int i = 0; i < Q; ++i) {
                    const int neighbour = neighbouring_fi_cell_index(x, y, z, i, length) * Q;
                    assert(neighbour % Q == 0); // Make sure it points to the start of a cell.
                    streamField[field_index + i] = collideField[neighbour + i];

                    assert(streamField[field_index + i] >= 0.0);
                }
            }
        }
    }
}