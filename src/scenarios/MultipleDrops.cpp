#include "MultipleDrops.hpp"
#include <LBMHelper.hpp>
#include "LBDefinitions.hpp"

MultipleDrops::MultipleDrops(ConfigParser &config){  }

void MultipleDrops::getFlagField(std::vector<flag_t> &flags, const coord_t &length) {

    const coord_t center_a= {5, 5, 30};
    const int radius_a = 3;

    const coord_t center_b = {25, 25, 20};
    const int radius_b = 8;

    const coord_t center_c = {25, 25, 36};
    const int radius_c = 3;

    const coord_t center_d = {40, 40, 20};
    const int radius_d = 8;

    const coord_t center_e = {8, 42, 33};
    const int radius_e = 6;

    const coord_t center_f = {40, 15, 25};
    const int radius_f = 7;

    const coord_t center_g = {8, 20, 15};
    const int radius_g = 6;

    const coord_t center_h= {30, 7, 35};
    const int radius_h = 4;

    const int waterHeight = 5;
    for (int z = 1; z < length[2] + 1; ++z) {
        for (int y = 1; y < length[1] + 1; ++y) {
            for (int x = 1; x < length[0] + 1; ++x) {
                flag_t entry = flag_t::EMPTY;
                double d_a = distance(x, y, z, center_a);
                double d_b = distance(x, y, z, center_b);
                double d_c = distance(x, y, z, center_c);
                double d_d = distance(x, y, z, center_d);
                double d_e = distance(x, y, z, center_e);
                double d_f = distance(x, y, z, center_f);
                double d_g = distance(x, y, z, center_g);
                double d_h = distance(x, y, z, center_h);

                if (d_a <= radius_a) {
                    entry = flag_t::FLUID;
                }
                if (d_b <= radius_b) {
                    entry = flag_t::FLUID;
                }
                if (d_c <= radius_c) {
                    entry = flag_t::FLUID;
                }
                if (d_d <= radius_d) {
                    entry = flag_t::FLUID;
                }
                if (d_e <= radius_e) {
                    entry = flag_t::FLUID;
                }
                if (d_f <= radius_f) {
                    entry = flag_t::FLUID;
                }
                if (d_g <= radius_g) {
                    entry = flag_t::FLUID;
                }
                if (d_h <= radius_h) {
                    entry = flag_t::FLUID;
                }

                if (z <= waterHeight) {
                    entry = flag_t::FLUID;
                }

                flags[indexForCell(x, y, z, length)] = entry;
            }
        }
    }

}

double MultipleDrops::distance(int x, int y, int z, const coord_t &center) {
    return std::sqrt(std::pow(center[0] - x, 2) + std::pow(center[1] - y, 2) +
                     std::pow(center[2] - z, 2));
}