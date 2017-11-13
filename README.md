# TSAR
## Time-Series computation Automation for Raster images
These scripts offer functionality to:
* compute multitemporal statistics
* apply transforming computations (e.g. time series spline smoothing)

over the third (i.e. time) dimension of 3D raster stacks.

Once loaded into memory, the 3D array is split along the y-axis into n parts of approximately equal size with n being the number of CPUs defined.
The parts are then passed to the CPUs for computation. Once finished the computation results are combined along the y-axis.

In case of smoothing the resulting array will have the same size as the original file.
In case of statistics computation the resulting array will be 2D or 3D with the third dimension representing the different computed statistics.

As a last step after writing the file the ENVI header file will be edited so that the band names are either the names of the original file (smoothing) or the names of the defined statistics functions (referred to as workers).

Input and Output format is a 3D ENVI binary stack.
