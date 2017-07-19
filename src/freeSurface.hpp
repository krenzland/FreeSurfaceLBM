#ifndef CFD_LAB_FREESURFACE_HPP
#define CFD_LAB_FREESURFACE_HPP

#include "LBDefinitions.hpp"
#include "LBMHelper.hpp"
#include "computeCellValues.hpp"
#include "streaming.hpp"
#include <assert.h>
#include <unordered_set>
#include <vector>
using gridSet_t = std::unordered_set<coord_t>;

std::array<double, 3> computeSurfaceNormal(const std::vector<double> &distributions,
                                           const std::vector<double> &mass, const coord_t &position,
                                           const coord_t &length);

void streamMass(const std::vector<double> &distributions, const std::vector<flag_t> &flags,
                std::vector<double> mass, const coord_t &length);

// Corresponds to section 4.3.
// TODO: Find better name!
void getPotentialUpdates(const coord_t &coord, double mass, double density, gridSet_t &filled,
                         gridSet_t &emptied);

void flagReinit(std::vector<double> distributions, std::vector<double> &mass,
                std::vector<flag_t> &flags, gridSet_t &filled, gridSet_t &emptied,
                const coord_t &length);

void distributeMass(const std::vector<double> &distributions, std::vector<double> &mass,
                    std::vector<flag_t> &flags, gridSet_t &filled, gridSet_t &emptied,
                    const coord_t &length);
#endif // CFD_LAB_FREESURFACE_HPP
