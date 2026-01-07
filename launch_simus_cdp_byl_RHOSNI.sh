#!/bin/bash

# List of experiments
experiments=(
    #RHOSNI_CROCUS
    #RHOSNI_CROCUS_R21
    #RHOSNI_CROCUS_GW1
    #RHOSNI_CROCUS_GW2
    #RHOSNI_SNOWPACK
    RHOSNI_SnowTran
    RHOSNI_SnowTran_Arctic
)

# List of sites
sites=(
    "SnowMIP/cdp"
    "SnowArctic/byl"
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
        if [ "${site}" == "SnowArctic/byl" ]; then
            run_prefix="run_peat"
        else
            run_prefix="run"
        fi

        # Copy the job_option_file from ref run
        cp inputFiles/${site}/job_options_file_${run_prefix}_TCZE_REF_TSNB_REF.txt inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt

        # Modify path in job_options_file
        sed -i "s/TCZE_REF_TSNB_REF/${experiment}/g" inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt

        # Create output file directory
        mkdir -p outputFiles/${site}/${run_prefix}_${experiment}

        # Copy and restart files from REF run
        cp inputFiles/${site}/${site_base}_init_${run_prefix}_TCZE_REF_TSNB_REF.nc inputFiles/${site}/${site_base}_init_${run_prefix}_${experiment}.nc
        cp inputFiles/${site}/rsfile_${run_prefix}_TCZE_REF_TSNB_REF.nc inputFiles/${site}/rsfile_${run_prefix}_${experiment}.nc

        # Launch simulation
        apptainer exec CLASSICcontainer.simg bin/CLASSIC_serial_${experiment} inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt 0/0 |& tee outputFiles/${site}/${run_prefix}_${experiment}/output.txt
    done
done


