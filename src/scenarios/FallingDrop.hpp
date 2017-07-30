//
// Created by lukas on 29/07/17.
//

#ifndef CFD_LAB_FALLINGDROP_HPP
#define CFD_LAB_FALLINGDROP_HPP

#include "Scenario.hpp"
#include <ConfigParser.hpp>

class FallingDrop : public Scenario {
  public:
    explicit FallingDrop(ConfigParser &config);
    void getFlagField(std::vector<flag_t> &flags, const coord_t &length) override;

  private:
    int waterHeight;
    int dropHeight;
    int dropRadius;
};

#endif // CFD_LAB_FALLINGDROP_HPP
