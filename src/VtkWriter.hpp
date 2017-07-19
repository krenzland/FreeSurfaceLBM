#ifndef CFD_LAB_VTKWRITER_HPP
#define CFD_LAB_VTKWRITER_HPP

#include "LBDefinitions.hpp"
#include "computeCellValues.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

class VtkWriter {
  public:
    VtkWriter(std::string filenameRoot, coord_t length, coord_t realLength, coord_t offset,
              coord_t processId);

    void write(const std::vector<double> &collideField, const std::vector<double> &mass,
               const std::vector<flag_t> &flagField, int t);

  private:
    std::string filenameRoot;
    coord_t length;
    coord_t realLength;
    coord_t offset;
    coord_t processId;
};

#endif // CFD_LAB_VTKWRITER_HPP
