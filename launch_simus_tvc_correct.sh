#!/bin/bash

# List of experiments
experiments=(
    PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne
)

# List of sites
sites=(
    "SnowArctic/tvc"
)

# Loop over each site
for site in "${sites[@]}"; do
    # Loop over each experiment in the list
    for experiment in "${experiments[@]}"; do
        # Print the current site and experiment
        echo "Processing site: ${site}, experiment: ${experiment}"
        
        # Extract the base name from the site path (e.g., "cdp" or "byl")
        site_base=$(basename ${site})

        # Determine the file name prefix based on site
        run_prefix="run_1peat"

        # Copy the job_option_file from ref run
        cp inputFiles/${site}/job_options_file_${run_prefix}_30min_ext.txt inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt

        # Modify path in job_options_file
        sed -i "s/1peat_30min_ext/1peat_${experiment}/g" inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt

        # Create output file directory
        mkdir -p outputFiles/${site}/${run_prefix}_${experiment}

        # Copy and restart files from REF run
        cp inputFiles/${site}/${site_base}_init_${run_prefix}_30min_ext.nc inputFiles/${site}/${site_base}_init_${run_prefix}_${experiment}.nc
        cp inputFiles/${site}/rsfile_${run_prefix}_30min_ext.nc inputFiles/${site}/rsfile_${run_prefix}_${experiment}.nc

        # Launch simulation
        apptainer exec CLASSICcontainer.simg bin/CLASSIC_serial_${experiment} inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt 0/0 |& tee outputFiles/${site}/${run_prefix}_${experiment}/output.txt
    done
done


