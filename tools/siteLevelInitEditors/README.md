Tools to edit exising site-level initialization
============

These tools are designed to make conversion or editing of site-level init files easier.

- `nc_to_cdl.py`: This script takes in a site-level init file and converts it the [`cdl` format](https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_utilities_guide.html) using `ncdump` tool (requires nco). It formats the cdl output in a manner that is easier to edit than the default ncdump format. After editing the resulting cdl a new netcdf init file can be easily created using `ncgen {file}.cdl -o {outfilename}.nc`.
- `convertInitToV1.1.py`: This script adds in the new variables needed to run CLASSIC v.1.1.x from files of the CLASSIC v1.0.x format.
- `expandv1.1Filefrom9PFTto12PFT.py`: Expands a CLASSIC v.1.1.x format file from 9 PFT to a 12 PFT. The same arrays are simply expanded allowing the use of the CLASSIC run parameters file for the 12 PFTs (which includes 2 shrub PFTs and sedges).
