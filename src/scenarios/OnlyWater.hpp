//
// Created by lukas on 29/07/17.
//

#ifndef CFD_LAB_ONLYWATER_HPP
#define CFD_LAB_ONLYWATER_HPP


#include "Scenario.hpp"

class OnlyWater : public Scenario {
public:
    void getFlagField(std::vector<flag_t> &flags, const coord_t &length) override;
};


#endif //CFD_LAB_ONLYWATER_HPP
