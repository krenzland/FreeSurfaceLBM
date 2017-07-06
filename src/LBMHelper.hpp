//
// Created by Andreas Reiser on 24.05.17.
//

#ifndef SOURCE_LBMHELPER_H
#define SOURCE_LBMHELPER_H

#include "LBDefinitions.hpp"

int indexForCell(int x, int y, int z, int *length);
int indexForCell(int x, int y, int z, const coord_t &length);
int inverseVelocityIndex(int index);

#endif // SOURCE_LBMHELPER_H
