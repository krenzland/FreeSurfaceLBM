#ifndef CFD_LAB_DAMBREAK_HPP
#define CFD_LAB_DAMBREAK_HPP

#include "Scenario.hpp"
#include <ConfigParser.hpp>

class DamBreak : public Scenario {
  public:
    explicit DamBreak(ConfigParser &config);
    void getFlagField(std::vector<flag_t> &flags, const coord_t &length) override;

  private:
    int damSize;
};

#endif // CFD_LAB_DAMBREAK_HPP
