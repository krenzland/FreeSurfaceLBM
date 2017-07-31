#include "VtkWriter.hpp"
#include <assert.h>

VtkWriter::VtkWriter(const std::string &filenameRoot, const coord_t &length)
    : filenameRoot(filenameRoot), length(length) {}

void VtkWriter::write(const std::vector<double> &collideField, const std::vector<double> &mass,
                      const std::vector<double> &density, const std::vector<flag_t> &flagField,
                      double stepSize, int t) {
    std::stringstream filenameBuilder;
    // Filename is of the format name_i_j_k.t, where ijk are the mpi
    filenameBuilder << filenameRoot << "." << t << ".vtk";
    auto filename = filenameBuilder.str();

    std::ofstream os(filename);
    if (!os) {
        throw std::runtime_error("Couldn't open file " + filename + ".");
    }

    // First write the header.
    os << "# vtk DataFile Version 2.0 \n"
       << "LBM data \n"
       << "ASCII \n\n"
       << "DATASET STRUCTURED_POINTS\n" // the grid type
       // Only the length without ghost cells is relevant here.
       << "DIMENSIONS " << length[0] + 2 << ' ' << length[1] + 2 << ' ' << length[2] + 2 << '\n'
       // Each process has a different origin.
       << "ORIGIN 0 0 0" << '\n'
       << "SPACING 1 1 1 \n\n";

    const auto noCells = (length[0] + 2) * (length[1] + 2) * (length[2] + 2);

    os << "POINT_DATA " << noCells << "\nVECTORS velocity float\n\n";

    double velocity[3] = {0};
    for (size_t i = 0; i < density.size(); ++i) {
        computeVelocity(&collideField[i * Q], density[i], velocity);
        os << velocity[0] / stepSize << ' ' << velocity[1] / stepSize << ' '
           << velocity[2] / stepSize << '\n';
    }

    os << "\nSCALARS density float 1"
       << "\nLOOKUP_TABLE default\n\n";

    for (size_t i = 0; i < density.size(); ++i) {
        os << density[i] << '\n';
    }

    os << "\nSCALARS cell_type int 1"
       << "\nLOOKUP_TABLE default\n\n";

    for (size_t i = 0; i < density.size(); ++i) {
        os << static_cast<int>(flagField[i]) << '\n';
    }

    os << "\nSCALARS mass float 1"
       << "\nLOOKUP_TABLE default\n\n";

    for (size_t i = 0; i < mass.size(); ++i) {
        os << mass[i] << '\n';
    }

    os << "\nSCALARS volumeOfFluid float 1"
       << "\nLOOKUP_TABLE default\n\n";

    for (size_t i = 0; i < mass.size(); ++i) {
        os << mass[i] / density[i] << '\n';
    }
}

void VtkWriter::writeMass(const std::vector<double> &mass, int t) {
    std::stringstream filenameBuilder;
    // Filename is of the format name_i_j_k.t, where ijk are the mpi
    filenameBuilder << filenameRoot << "." << t << ".vtk";
    auto filename = filenameBuilder.str();

    std::ofstream os(filename);
    if (!os) {
        throw std::runtime_error("Couldn't open file " + filename + ".");
    }

    // First write the header.
    os << "# vtk DataFile Version 2.0 \n"
       << "LBM data \n"
       << "ASCII \n\n"
       << "DATASET STRUCTURED_POINTS\n" // the grid type
       // Only the length without ghost cells is relevant here.
       << "DIMENSIONS " << length[0] + 2 << ' ' << length[1] + 2 << ' ' << length[2] + 2 << '\n'
       // Each process has a different origin.
       << "ORIGIN 0 0 0" << '\n'
       << "SPACING 1 1 1 \n\n";

    os << "\nSCALARS mass float 1"
       << "\nLOOKUP_TABLE default\n\n";

    for (size_t i = 0; i < mass.size(); ++i) {
        os << mass[i] << '\n';
    }
}
