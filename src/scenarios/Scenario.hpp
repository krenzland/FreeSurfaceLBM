#ifndef CFD_LAB_SCENARIO_HPP
#define CFD_LAB_SCENARIO_HPP

#include <LBDefinitions.hpp>
#include <vector>

class Scenario {
public:
    // Initialises the flag field without the boundary layer.
    virtual void getFlagField(std::vector<flag_t> &flags, const coord_t &length)= 0;
};


#endif //CFD_LAB_SCENARIO_HPP
