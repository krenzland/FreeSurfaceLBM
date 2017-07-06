#ifndef _STREAMING_H_
#define _STREAMING_H_

#include "LBDefinitions.hpp"
#include <vector>

/** carries out the streaming step and writes the respective distribution
 * functions from collideField to streamField.
 */
void doStreaming(const std::vector<double> &collideField, std::vector<double> &streamField,
                 const std::vector<flag_t> &flagField, const coord_t &length);
#endif
