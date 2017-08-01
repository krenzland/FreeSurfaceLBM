# Compilation

Compile using `cmake . ` followed by `make`.
Cmake should find OpenMP automatically, at least for GCC.

# Usage
Use one of the supplied parameter files.
You can find a selection of nice examples in the examples/ folder.

Note that we use an OpenMP based.
This means that the number of used threads is read from the environment.
To set it, simply run `export OMP_NUM_THREADS=4`.
The program is used by simply running `lbsim params.dat`, where params.dat is a parameter file.
The results can then be found in the folder results/.

# Rendering
There are two ways the results can be rendered in ParaView.
- Use the filter contour, and set the contour value to 0.5 of the volumeOfFluid variable.
- Use volume rendering and colour it with the volumeOfFluid variable.

# Documentation
The report is in the folder docs/.