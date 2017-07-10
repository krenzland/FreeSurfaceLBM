#ifndef CFD_LAB_FREESURFACE_HPP
#define CFD_LAB_FREESURFACE_HPP

#include <vector>
#include <assert.h>
#include <unordered_set>
#include "computeCellValues.hpp"
#include "LBDefinitions.hpp"
#include "LBMHelper.hpp"
#include "streaming.hpp"
using gridSet_t = std::unordered_set<coord_t>;

void streamMass(const std::vector<double> &distributions, const std::vector<flag_t> &flags, std::vector<double> mass,
                const coord_t &length);

// TODO: Boundary Conditions for Interface cells, maybe perform during streaming step?

// Corresponds to section 4.3.
// TODO: Find better name!
void getPotentialUpdates(const coord_t &coord, double mass, double density, gridSet_t &filled, gridSet_t &emptied);

void flagReinit(std::vector<double> &mass, std::vector<flag_t> &flags, gridSet_t &filled, gridSet_t &emptied,
                const coord_t &length);

#endif //CFD_LAB_FREESURFACE_HPP
