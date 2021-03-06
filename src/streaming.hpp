#ifndef _STREAMING_H_
#define _STREAMING_H_

#include "LBDefinitions.hpp"
#include <vector>

/** carries out the streaming step and writes the respective distribution
 * functions from collideField to streamField.
 */
void doStreaming(const std::vector<double> &collideField, std::vector<double> &streamField, const std::vector<double> &mass,
                 const coord_t &length, const std::vector<flag_t> &flagField, std::vector<neighborhood_t> &neighborhood,
                 std::vector<double> &fluidFraction);
int neighbouring_fi_cell_index(int x, int y, int z, int fi, const coord_t &length);
#endif
