# Creation and formatting of model input files {#makeInputFiles}

1. @ref makeMet
2. @ref makeInit
3. @ref ghgfiles
4. @ref makeOther
5. @ref inputFileForm

---

# Preparation of files for meteorological inputs {#makeMet}

CLASSIC requires the meteorological inputs as described here. The input files are netCDF format. A netcdf dump of the files will look like this:

        netcdf pres_v1.1.5_T63_chunked_1700_2017 {
        dimensions:
            lon = 128 ;
            lat = 64 ;
            time = UNLIMITED ; // (464280 currently)
        variables:
            float lon(lon) ;
                lon:standard_name = "longitude" ;
                lon:long_name = "longitude" ;
                lon:units = "degrees_east" ;
                lon:axis = "X" ;
                lon:_Storage = "contiguous" ;
                lon:_Endianness = "little" ;
            float lat(lat) ;
                lat:standard_name = "latitude" ;
                lat:long_name = "latitude" ;
                lat:units = "degrees_north" ;
                lat:axis = "Y" ;
                lat:_Storage = "contiguous" ;
                lat:_Endianness = "little" ;
            double time(time) ;
                time:standard_name = "time" ;
                time:units = "**day as %Y%m%d.%f**" ;
                time:calendar = "proleptic_gregorian" ;
                time:axis = "T" ;
                time:_Storage = "chunked" ;
                time:_ChunkSizes = 1 ;
                time:_Endianness = "little" ;
            float pres(time, lat, lon) ;
                pres:long_name = "Pressure" ;
                pres:units = "Pa" ;
                pres:_FillValue = 9.96921e+36f ;
                pres:missing_value = 9.96921e+36f ;
                pres:_Storage = "chunked" ;
                pres:_ChunkSizes = 464280, 8, 16 ;
                pres:_Endianness = "little" ;

        }
Importantly,

- Only one variable per file besides lon, lat, and time.
- Note the time units.
- The file is chunked (which depends on the grid being used. What you see here is optimal for T63 global runs). Chunking is not needed for site-level runs. More on chunking in @ref inputFileForm.
- CLASSIC expects the first time step to be 0 hour 0 minute (in hh:mm:ss format 00:00:00).

If you have existing ACSCII met files from a site, see [ASCII to NetCDF met file loader](@ref asciiMet) to use the provided tool to convert them to the appropriate netCDF format.

# Preparation of the model initialization file {#makeInit}

If you have the old format .INI (and CTEM's .CTM) initialization files, there is a [tool to convert them to netCDF format](@ref initTool) for use in CLASSIC (only up to version 1.0!). The tool itself is located in tools/initFileConverter.

If you have a netCDF format initialization/restart file that you wish to edit, you can use the [script created for that purpose](@ref modifyRS). It is located in tools/modifyRestartFile. This is designed for regional/global scale edits. It may be simpler to just create your own script using a Jupyter notebook with the xarray library, just depends on your comfort-level manipulating these files. For site-level files (single point) it is likely easier to use a combination of `ncdump`/`ncgen`. There are tools to assist this (especially tools/siteLevelInitEditors/nc_to_cdl.py)

To see properly formatted files, perform an `ncdump` on any of the files in the CBC. Those initialization files are setup for a biogeochemistry (CTEM on) run with peatlands, nitrogen cycle, and PFT competition variables included. Note that for a {physics-only run (CTEM off), no interactive N cycle, no peatlands, no competition}, many of the variables will not be read in/required, similarly for a CTEM on run some CLASS-only variables are not required. See the [manual's mainpage](@ref main) for links to sections describing the variables.

# Greenhouse gas inputs files {#ghgfiles}

If you have existing ACSCII GHG files, see [our tool to create GHG input files](@ref makeGHGfiles) to convert them to the appropriate netCDF format.

# Making input files for other input variables {#makeOther}

Most other CLASSIC input files are not desired, or needed, for point runs of the model. If you require input files for regional or global simulations we may be able to provide you with versions we use in our runs. The possible inputs are listed in @ref CTEMaddInputs depending on model configuration. Because these files are not required of most users we have not set up tools to help generate the files.

# Some notes on input file format {#inputFileForm}

These points below are most relevant for those running globally or at high resolution regionally. They will have minimal impact upon site-level runs. 

**Input File Format**

Because CLASSIC utilizes the parallel version of netCDF-4 libraries, all input files must be in netCDF-4 format. Many data files continue to be distributed in netCDF "classic" format, sometimes referred to as netCDF-3. There are various tools available to convert from one format to the other. One such utilitiy is `ncks` which is part of the [nco package](http://nco.sourceforge.net).  For example:
```
ncks -4 input_file output_file
```
Further information regarding netCDF can be found on Unidata's [FAQ](https://www.unidata.ucar.edu/software/netcdf/docs/faq.html). 

> **Simply converting to netCDF-4 format is not usually all that is required to obtain good I/O performance when running the parallel version of CLASSIC. Please note carefully the following section.**

**Chunking Input Files**

One of the requirements of the netCDF-4 format is that if a record dimension is present, the file must be chunked. So when a netCDF classic file with an unlimited time dimension is converted to netCDF-4 format, it is automatically "chunked". Chunking is meant to enhance I/O performance, but this is only the case if the chunks are chosen in accordance to how the file is to be accessed. See this link for details of why [chunking matters](https://www.unidata.ucar.edu/blogs/developer/en/entry/chunking_data_why_it_matters). `It is important to note that a poor choice of chunk sizes can, in fact, lead to very poor I/O performance. Popular utilities such as "ncks" and "cdo" typically do not, by default, chunk netCDF-4 files in a way that works well with CLASSIC.`

The parallel version of CLASSIC runs the entire simulation at each grid cell in a domain on an individual processor. This means that all time steps in the input file must be read in for that grid cell at the very start of the simulation. In order to ensure good read performance, the input files should be chunked so that all the time steps for a given grid cell are included in the same chunk. However, it is not recommended to chunk a file so that each chunk only contains the time steps at a single grid point as it would mean very poor read performance when accessing the file as 2D grids in latitude/longitude. The best compromise is to chunk the files so that each chunk contains all time steps for multiple grid points. For example, using ncks on a netCDF file of T63 global grids with a `time` dimension consisting of 39420 time steps, a `lat` dimension of 64 latitudes and a `lon` dimension of 128 longitudes, a suitable command line would be: 
```
ncks -h -4 --cnk_plc=g3d --cnk_dmn=time,39420 --cnk_dmn=lat,8 --cnk_dmn=lon,16 input_file output_file
```
And for the high resolution Canada domain (310x160 rotated lat-lon grid):
```
ncks -h -4 --cnk_plc=g3d --cnk_dmn=time,39420 --cnk_dmn=lat,10 --cnk_dmn=lon,10 input_file output_file
```
It is recommended that:

1.  All processing of the input file with other utilities such as `cdo` or other `nco` tools be completed before this final step since those utilities may rechunk the file using different chunk sizes.
2.  The input file is in netCDF classic format so that the chunking and conversion to netCDF-4 format (via the "-4" switch) is done at the same time.
3.  The entire number of time steps in the input file is used for the chunk size in the time dimension. 
4.  The chunk sizes for the lat and lon dimensions in the above examples be used for those domains. However, it may be that for files with longer time series, the chunk sizes for the lat and lon dimensions need to be adjusted (i.e. reduced).
5.  The output file is checked to ensure that the chunking was done properly. This can be done using `ncdump -hs output_file`. The output from `ncdump` should show the variable with the chunk sizes used by the `ncks` utility and the time variable should also be chunked with the total number of time steps. The lat and lon variables will not be chunked (i.e. will be shown as "contiguous" in the ncdump output).

> **Caveats for ECCC users**

> **Local implementations of nco/ncks may not support all the features needed to properly process and chunk files suitable for use on the ECCC Science network. On each ppp, /usr/bin/ncks appears to work well and it is recommended to use it rather than a local implementation of nco to do the chunking if there are issues encountered locally.**

> **Chunking (via ncks) requires that the entire file be placed into memory. If there is insufficient RAM, the chunking will fail. Note that there is about 187 GB of RAM available on each node of the ppp5/6 machines.**
