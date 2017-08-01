#ifndef CFD_LAB_MULTIPLEDROPS_HPP
#define CFD_LAB_MULTIPLEDROPS_HPP

#include "Scenario.hpp"
#include <ConfigParser.hpp>

class MultipleDrops : public Scenario {
public:
    explicit MultipleDrops(ConfigParser &config);
    void getFlagField(std::vector<flag_t> &flags, const coord_t &length) override;

private:
    double distance(int x, int y, int z, const coord_t &center);
};
#endif //CFD_LAB_MULTIPLEDROPS_HPP
