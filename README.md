# Scripts to process data collected with the IR screen setup
The setup is described in the article Haugeneder, M., Lehning, M., Reynolds, D. et al. A Novel Method to Quantify Near-Surface Boundary-Layer Dynamics at Ultra-High Spatio-Temporal Resolution. Boundary-Layer Meteorol 186, 177–197 (2023). https://doi.org/10.1007/s10546-022-00752-3.

This is a reorganization of the scripts

## Steps to process the raw *.irb-sequences to netCDF4-files

1. convert *.irb-files to *.txt-files (ASCII, one file per frame) with IRB2ASCII (separate software). Check that decimal seperator is "." and space between values is "tab". Multiple sequences can be converted to a single folder (just not e.g. two times "irdata00")
2. Use script code/process_irraw.jl to convert single .txt-files to one netCDF4-file per sequence. Change variable "pathtofile" to point to folder containing *.txt-files. Output is written to same folder.

## Raw data visualization

The raw data can be visualized using the code/check_view_data.jl script.
