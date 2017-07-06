#include "VtkWriter.hpp"

VtkWriter::VtkWriter(std::string filenameRoot, coord_t length, coord_t realLength, coord_t offset,
                     coord_t processId)
    : filenameRoot(filenameRoot), length(length), realLength(realLength), offset(offset),
      processId(processId) {}

void VtkWriter::write(const std::vector<double> &collideField, const std::vector<flag_t> &flagField,
                      int t) {
    std::stringstream filenameBuilder;
    // Filename is of the format name_i_j_k.t, where ijk are the mpi
    filenameBuilder << filenameRoot << '_' << processId[0] << '_' << processId[1] << '_'
                    << processId[2] << "_." << t << ".vtk";
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
       << "DIMENSIONS " << realLength[0] << ' ' << realLength[1] << ' ' << realLength[2] << '\n'
       // Each process has a different origin.
       << "ORIGIN " << offset[0] << ' ' << offset[1] << ' ' << offset[2] << '\n'
       << "SPACING 1 1 1 \n\n";

    // We precompute the density for every cell.
    const auto number_cells = (length[0] + 2) * (length[1] + 2) * (length[2] + 2);
    const auto numberRealCells = realLength[0] * realLength[1] * realLength[2];

    auto density = std::vector<double>(number_cells);
    for (size_t i = 0; i < density.size(); ++i) {
        density[i] = computeDensity(&collideField[i * Q]);
    }

    os << "POINT_DATA " << numberRealCells << "\nVECTORS velocity float\n\n";

    double velocity[3] = {0};
    for (size_t i = 0; i < density.size(); ++i) {
        if (flagField[i] != flag_t::PARALLEL_BOUNDARY) {
            computeVelocity(&collideField[i * Q], density[i], velocity);
            os << velocity[0] << ' ' << velocity[1] << ' ' << velocity[2] << '\n';
        }
    }

    os << "\nSCALARS density float 1"
       << "\nLOOKUP_TABLE default\n\n";

    for (size_t i = 0; i < density.size(); ++i) {
        if (flagField[i] != flag_t::PARALLEL_BOUNDARY) {
            os << density[i] << '\n';
        }
    }

    os << "\nSCALARS cell_type int 1"
       << "\nLOOKUP_TABLE default\n\n";

    for (size_t i = 0; i < density.size(); ++i) {
        if (flagField[i] != flag_t::PARALLEL_BOUNDARY) {
            os << static_cast<int>(flagField[i]) << '\n';
        }
    }
}
