# Contour Boxplot

Computing the contour boxplot for an ensmeble of 2D isocontours. This code has been tuned for [SREF ensemble dataset (text file format)](http://www.nco.ncep.noaa.gov/pmb/products/sref/). For more details regarding the contour boxplot which is a nonparamteric technique to summarize the main features of an ensmeble of isocontour consult the [contour boxplot paper](http://www.cs.miami.edu/home/mirzargar/papers/contour_boxplot.pdf) or its [wiki page](https://en.wikipedia.org/wiki/Contour_boxplot).

## Requirements
```
CMake 2.8+
C++11 compiler
ITK4.8+
```
## Building
```
cd contour_boxplot
mkdir build
cd build
cmake ../src/ or cmake ../src/ -DITK_DIR=<ITK build directory> 
make
```

## Running
* Note that this code is the original version (not necessary fast or memory optimized). The SVG generator is also tuned for [SREF ensemble](http://www.nco.ncep.noaa.gov/pmb/products/sref/).

If you run the code with no parameters there are more details about assumptions in this code that might need to be changed.

```
./contourboxplot <src_dir> <ens_size> <level_set> <des_dir> <Xdim> <Ydim>
```

Note:
* files are considered to be linearized - if change needed consult function ```LoadFieldFileWithoutCoord``` in EnsembleHandler.cpp
* the files should be named as mem_c_0 to mem_c_ens_size - if change needed consult ```LoadFieldEnsemble``` in EnsembleHandler.cpp
* the value of j has been hard coded in the main file (```j=2```) - it also has been hard coded partially in ```getTimePercents``` function in CBD.cpp
* the outlier num is specified as the number of members with zero depth value. If every member has a non-zero depth value then it has been 
