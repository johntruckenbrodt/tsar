# TSAR
## Time-Series computation Automation for Raster data
These scripts offer functionality to:
* compute multitemporal statistics
* apply transforming computations (e.g. time series spline smoothing)

over the third (i.e. time) dimension of 3D raster stacks.

Input is a raster stack which might be pointing to a single file or virtually to a list fo files with same spatial extent.
A master process divides the number of lines by the number of parallel processes and passes each instructions on which subset to read from the stack. Each process then reads its subset and computes a set of functions (workers) over the time dimension.

Once finished the computation results are gathered by the master process and combined along the y-axis.

The defined workers take a vector as input (i.e. a pixel time series) and can output a single value (e.g. a multi-temporal statistic) or several values (e.g. time series smoothing).
The number of resulting bands is checked for a single pixel before the parallel processing is started to determine whether any of the resulting files might already exist and deletes single worker functions from the list if existing files are not to be overwritten.

The names of the computed bands can either be explicitely defined by the user or are otherwise generated from the names of the workers.
In this case, if individual workers return several values, each band will carry numeric suffix (e.g. worker1_1, worker1_2, worker2_1, etc.).

Depending on whether the output is to be a single file or separate files, an ENVI stack or several GeoTiff files are written.
In case an ENVI stack is written, the file itself is named as defined by out.name and the bandnames are written to the HDR file.
If individual GeoTiffs are written, out.name is going to be a directory with the bandnames being the names of the GeoTiff files.
