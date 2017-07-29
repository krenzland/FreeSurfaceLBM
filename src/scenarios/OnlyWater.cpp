#include "OnlyWater.hpp"

void OnlyWater::getFlagField(std::vector<flag_t> &flags, const coord_t &length) {
    for (auto &flag : flags) {
        flag = flag_t::FLUID;
    }
}
