# Scripts to Launch CLASSIC on the General Purpose Science Cluster for Collaboration (gpscc3).

The following scripts may be used to compile, run and monitor CLASSIC.

* make_classic_collab.sh
* classic_submit_collab.sh
* jobstat_collab

The following script is normally run automatically to stitch output files from multiple nodes onto the simulation domain.
Users shouldn't have to run it (but could do so if necessary/desired). 

* classic_stitch_netcdf.sh

It is included here:

[https://gitlab.science.gc.ca/CLASSIC_development/generalTools/tree/master/classic_run_scripts](https://gitlab.science.gc.ca/CLASSIC_development/generalTools/tree/master/classic_run_scripts)

and works on `gpscc3` as well as on the HPC clusters.

## Compiling CLASSIC

It is strongly recommended that `make_classic_collab.sh` be used when compiling CLASSIC for `gpscc3`.
The script sets required variables and library paths prior to running `make` on the `Makefile`.
The following commands will compile an executable for CLASSIC:

``` 
cd <directory containing the Makefile for CLASSIC>
<path to the script>/make_classic_collab.sh
```

A copy of `make_classic_collab.sh` is maintained in the directory: `~rec001/public/classic`.
If the script is used often, a symbolic link can be created in the user's own $HOME/bin directory, e.g.:

```
cd ~/bin
ln -s ~rec001/public/classic/make_classic_collab.sh make_classic.sh 

```

## Submitting CLASSIC

1. Create an output directory for the run (under `/space/hall0/work/eccc/crd/ccrp/$USER`) and place 
   at least an edited copy of the job options file into it. 
   It is usually recommended to also place a copy of the run parameters and variable descriptor files used for the run
   (specifying the path to these in the job options file). It is also recommended to include a copy of a restart 
   (or initialization file), again with the path specified in the job options file (i.e. `rs_file_to_overwrite`).
   Note that, contrary to the description in the job options file, the file is not actually overwritten.
   Rather, a copy of it is made on each node where the model is run and those copies are overwritten.

2. Obtain a copy of `classic_submit_collab.sh` and ensure that the following parameters are set correctly:
```
   output_directory=<full path to the run directory>
   job_options_file=<full path to the job options file>
   executable=<full path to the compiled CLASSIC executable>
   email=<user's email address to send job status notifications>
```
3. Launch the run by executing the script (without parameters) and follow the instructions on the screen.
   After the job is finished, the stitched model output netCDF files will be located under the `netcdf_files`
   directory. The temporary `netcdf_files_XX` (XX=01-99) directories may be safely deleted once the run is confirmed to
   have completed successfully.

## Monitoring the Job

The status of the run can be checked by running `jobstat_collab` on the system where the job was submitted. Optional parameters
allow jobs to checked on other systems and/or for other users.

A copy is also maintained in the directory: `~rec001/public/classic`.
If the script is used often, a symbolic link can be created in a user's own $HOME/bin directory, e.g.:

```
cd ~/bin
ln -s ~rec001/public/classic/jobstat_collab jobstat
```

## AMBER and Edifice

AMBER and Ediface are not supported yet on `gpscc3`. This section will be updated if/when these are made available.
