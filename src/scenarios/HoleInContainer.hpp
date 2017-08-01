#ifndef CFD_LAB_HOLEINCONTAINER_HPP
#define CFD_LAB_HOLEINCONTAINER_HPP

#include "Scenario.hpp"
#include <ConfigParser.hpp>

class HoleInContainer : public Scenario {
public:
    explicit HoleInContainer(ConfigParser &config);
    void getFlagField(std::vector<flag_t> &flags, const coord_t &length) override;
};
#endif //CFD_LAB_HOLEINCONTAINER_HPP
