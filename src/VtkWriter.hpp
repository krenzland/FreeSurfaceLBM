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
    VtkWriter(const std::string &filenameRoot, const coord_t &length);
    void write(const std::vector<double> &collideField, const std::vector<double> &mass,
                   const std::vector<flag_t> &flagField, double stepSize, int t, std::vector<double> &fluidFraction);
    void writeMass(const std::vector<double> &mass, int t);

  private:
    std::string filenameRoot;
    coord_t length;
};

#endif // CFD_LAB_VTKWRITER_HPP
