
#ifndef CFD_LAB_MIDDLEWALL_HPP
#define CFD_LAB_MIDDLEWALL_HPP

#include "Scenario.hpp"
#include <ConfigParser.hpp>

class MiddleWall : public Scenario {
public:
    explicit MiddleWall(ConfigParser &config);
    void getFlagField(std::vector<flag_t> &flags, const coord_t &length) override;
private:
    int wallLength;
};

#endif //CFD_LAB_MIDDLEWALL_HPP
