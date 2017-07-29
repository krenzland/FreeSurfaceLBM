#include <LBMHelper.hpp>
#include "FallingDrop.hpp"

FallingDrop::FallingDrop(ConfigParser &config) {
    dropHeight = config.parse<int>("dropHeight");
    dropRadius = config.parse<int>("dropRadius");
    waterHeight = config.parse<int>("dropWaterHeight");
}

void FallingDrop::getFlagField(std::vector<flag_t> &flags, const coord_t &length) {
    const int middleX = length[0] / 2;
    const int middleY = length[1] / 2;
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = 1; x < length[0] + 1; ++x) {
                flag_t entry = flag_t::EMPTY;
                const double distance = std::sqrt(
                        std::pow(middleX - x, 2) + std::pow(middleY - y, 2) + std::pow(dropHeight - z, 2));
                if (distance <= dropRadius) {
                    entry = flag_t::FLUID;
                }
                if (z <= length[2]/10) {
                    entry = flag_t::FLUID;
                }
                flags[indexForCell(x, y, z, length)] = entry;
            }
        }
    }
}

