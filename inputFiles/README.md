# Benchmarking with FLUXNET sites

## Philosophy

FLUXNET sites present an excellent source of observation-based data to evaluate CLASSIC. Ideally evaluation occurs at numerous sites to moderate the impact of site-level characteristics that can make it difficult for the model to reproduce observations. To facilitate this CLASSIC has a benchmarking suite that drives the model with the meteorological observations from the site and then evaluates the model outputs with the tower estimates. Over time the number of site is expected to increase and also to be refined. This suite has the advantage of each site being explicitly setup to try and best emulate conditions at each site and the model is driven with tower-observed meteorology.

## Suggested workflow

We use git to keep track of the model initialization files for each site to ensure we are able to match model input files to the code version. It is also important to keep track of changes in the model initialization files that may impact model results so these changes are not attributed, incorrectly, to changes in the model code. 

Since git is not well-suited to track binary formats we convert each model initialization file to a text format (cdl) using nco tools (available in the CLASSIC Singularity container). We also retain a yaml file for each site for model setup and notes about the site. Both the cdl and yaml files are tracked in our git repository. To facilitate working with the cdl and netcdf file formats we have the following tools:
- `createNetcdfsFromCDLs.py`: After cloning the repository this script can be run to convert all of the sites' cdl files to the netcdf format, which is needed to run the model at each site. 
- `updateRepository.py`: This converts the model netcdf initialization files to cdl format and can also do a git commit to the repository.

The other files necessary to run the model (CO2, meteorological, N-related) can be found either via the CLASSIC Community Zenodo page or in the scrd530 account (internal to ECCC only). These files are not tracked via this git repository due to their size.

## Adding new tower sites

If you are adding new sites to the benchmarking suite, the easiest solution is to run `/tools/siteLevelInitEditors/nc_to_cdl.py` on an existing netcdf. Then create a new folder for the site, copy that newly created cdl file and a yaml file (renaming appropriately) into the new folder. Edit the cdl and yaml files to be correct for the new site then run `ncgen {infile}.cdl -o {outfile}.nc` to generate a netcdf for model use. Make sure also to make the meteorology available to other model users.