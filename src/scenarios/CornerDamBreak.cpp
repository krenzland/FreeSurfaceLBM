#include "CornerDamBreak.hpp"
#include <LBMHelper.hpp>

CornerDamBreak::CornerDamBreak(ConfigParser &config) { damSize = config.parse<int>("damSize"); }

void CornerDamBreak::getFlagField(std::vector<flag_t> &flags, const coord_t &length) {
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = 1; x < length[0] + 1; ++x) {
                if (x< damSize && y < damSize) {}
                else {
                    flags[indexForCell(x, y, z, length)] = flag_t::EMPTY;
                }

            }
        }
    }
}