#include "HoleInContainer.hpp"
#include <LBMHelper.hpp>

HoleInContainer::HoleInContainer(ConfigParser &config) {  }

void HoleInContainer::getFlagField(std::vector<flag_t> &flags, const coord_t &length) {
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = 1; x < length[0] + 1; ++x) {
                if (z >= length[2]/2+2 && x <= length[0]/2) {}
                else {
                    flags[indexForCell(x, y, z, length)] = flag_t::EMPTY;
                }


            }
        }
    }

    //vertical wall
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            if (y >= length[1]/2-2 && y <= length[1]/2+2 &&
                z >= length[2]/2+2 && z <= length[2]/2+6) {}
                //hole 5X5
            else {
                flags[indexForCell(length[0]/2, y, z, length)] = flag_t::NO_SLIP;
                flags[indexForCell(1+length[0]/2, y, z, length)] = flag_t::NO_SLIP;
            }

        }
    }

    //horizontal wall
    for (int y = 1; y < length[1] + 1; ++y) {
        for (int x = 1; x < length[0] / 2 ; ++x) {
            flags[indexForCell(x, y, length[2]/2, length)] = flag_t::NO_SLIP;
            flags[indexForCell(x, y, 1+length[2]/2, length)] = flag_t::NO_SLIP;
        }
    }
}