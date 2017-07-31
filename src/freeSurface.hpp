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

std::array<double, 3> computeSurfaceNormal(const std::vector<double> &fluidFraction, const coord_t &length, const std::vector<double> &mass,
                                           const coord_t &position, const std::vector<flag_t> &flags);

void streamMass(const std::vector<double> &distributions, std::vector<double> &fluidFraction, const coord_t &length,
                std::vector<double> &mass, const std::vector<neighborhood_t> &neighborhood,
                const std::vector<flag_t> &flags);

double calculateSE(const std::vector<double> &distributions, const std::vector<flag_t> &flags,
                   const coord_t &coord, const coord_t &length, const int curFiIndex,
                   const std::vector<neighborhood_t> &neighborhood);
// Corresponds to section 4.3.
// TODO: Find better name!
void getPotentialUpdates(const std::vector<double> &mass, std::vector<double> &fluidFraction, std::vector<flag_t> &flags,
                         std::vector<neighborhood_t> &neighborhood, const coord_t &length);
void flagReinit(std::vector<double> &distributions, std::vector<double> &mass, std::vector<double> &fluidFraction,
                const coord_t &length, std::vector<flag_t> &flags);

void distributeMass(const std::vector<double> &distributions, std::vector<double> &mass, const coord_t &length,
                    std::vector<flag_t> &flags, std::vector<double> &fluidFraction);
#endif // CFD_LAB_FREESURFACE_HPP
