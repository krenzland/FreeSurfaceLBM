#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

#include "LBDefinitions.hpp"
#include <vector>

/** handles the boundaries in our simulation setup */
void treatBoundary(std::vector<double> &collideField, const std::vector<flag_t> &flagField,
                   const boundary_t &boundaryConditions, const coord_t &length);

#endif
