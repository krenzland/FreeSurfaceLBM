#include "DamBreak.hpp"
#include <LBMHelper.hpp>

DamBreak::DamBreak(ConfigParser &config) { damSize = config.parse<int>("damSize"); }

void DamBreak::getFlagField(std::vector<flag_t> &flags, const coord_t &length) {
    // Breaking  dam
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = damSize; x < length[0] + 1; ++x) {
                flags[indexForCell(x, y, z, length)] = flag_t::EMPTY;
            }
        }
    }
    for (int y = 1; y < length[1] + 1; ++y) {
        for (int x = 1; x < length[0] + 1; ++x) {
            flags[indexForCell(x, y, length[2], length)] = flag_t::EMPTY;
        }
    }

//    for (int z = 1; z < length[1] + 1; ++z) {
//        for (int x = 1; x < length[0] + 1; ++x) {
//            flags[indexForCell(x, 1, z, length)] = flag_t::EMPTY;
//            flags[indexForCell(x, length[1], z, length)] = flag_t::EMPTY;
//        }
//    }


    for (int z = 1; z < length[2]/3; ++z) {
        for (int y = 1; y < length[1]/2; ++y) {
            flags[indexForCell(length[0]/2, y, z, length)] = flag_t::NO_SLIP;
            flags[indexForCell(1+length[0]/2, y, z, length)] = flag_t::NO_SLIP;
        }
    }
}
