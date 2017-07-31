#ifndef CFD_LAB_TIMESTEP_HPP
#define CFD_LAB_TIMESTEP_HPP

#include "LBDefinitions.hpp"
#include "computeCellValues.hpp"
#include <array>
#include <vector>

// TODO: Use openMP!
// TODO: Maybe modify inflow velocity!
// Returns pair of tau and time step.
std::pair<double, double> adaptTimestep(std::vector<double> &distributions,
                                        std::vector<double> &density, std::vector<double> &mass,
                                        std::vector<flag_t> &flags,
                                        std::array<double, 3> &gravitation, double oldTimeStep,
                                        double oldTau, double smagorinskyConstant,
                                        bool allowIncrease);

#endif // CFD_LAB_TIMESTEP_HPP
