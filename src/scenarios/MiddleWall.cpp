//
// Created by Andreas Reiser on 31.07.17.
//

#include "MiddleWall.hpp"
#include <LBMHelper.hpp>

MiddleWall::MiddleWall(ConfigParser &config) { wallLength = config.parse<int>("wallLength"); }

void MiddleWall::getFlagField(std::vector<flag_t> &flags, const coord_t &length) {
    int yCoordWall = length[1] / 2;

    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = 1; x < length[0] + 1; ++x) {
                if (x < (length[0]/2)-2 && y < yCoordWall) {}
                else {
                    flags[indexForCell(x, y, z, length)] = flag_t::EMPTY;
                }

            }
        }
    }

    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = yCoordWall; y < yCoordWall+2; ++y) {
            for (int x = 0; x < wallLength+1; ++x) {
                flags[indexForCell(x, y, z, length)] = flag_t::NO_SLIP;
            }
        }
    }
}