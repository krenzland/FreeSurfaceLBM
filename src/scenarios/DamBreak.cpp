#include <LBMHelper.hpp>
#include "DamBreak.hpp"

DamBreak::DamBreak(ConfigParser &config) {
    damSize = config.parse<int>("damSize");
}

void DamBreak::getFlagField(std::vector<flag_t> &flags, const coord_t &length) {
    // Breaking  dam
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = damSize; x < length[0] + 1; ++x) {
                flags[indexForCell(x, y, z, length)] = flag_t::EMPTY;
            }
        }
    }
}

