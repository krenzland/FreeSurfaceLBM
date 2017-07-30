#ifndef CFD_LAB_CORNERDAMBREAK_H
#define CFD_LAB_CORNERDAMBREAK_H

#include "Scenario.hpp"
#include <ConfigParser.hpp>

class CornerDamBreak : public Scenario {
public:
    explicit CornerDamBreak(ConfigParser &config);
    void getFlagField(std::vector<flag_t> &flags, const coord_t &length) override;

private:
    int damSize;
};

#endif //CFD_LAB_CORNERDAMBREAK_H
