#!/bin/bash

# List of experiments
list=(
    WSNCAP_0
    WSNCAP_20
)

# Loop over each experiment in the list
for experiment in "${list[@]}"; do
    # Print the current experiment
    echo "Processing experiment: ${experiment}"
    
    # Copy the job_option_file from ref run
    cp inputFiles/SnowMIP/cdp/job_options_file_run_TCZE_REF_TSNB_REF.txt inputFiles/SnowMIP/cdp/job_options_file_run_${experiment}.txt

    # Modify path in job_options_file
    sed -i "s/TCZE_REF_TSNB_REF/${experiment}/g" inputFiles/SnowMIP/cdp/job_options_file_run_${experiment}.txt
    sed -i "s/default_run_parameters_shrubs/default_run_parameters_shrubs_${experiment}/g" inputFiles/SnowMIP/cdp/job_options_file_run_${experiment}.txt

    # Create output file directory
    mkdir -p outputFiles/SnowMIP/cdp/run_${experiment}

    # Copy and restart files from REF run
    cp inputFiles/SnowMIP/cdp/cdp_init_run_TCZE_REF_TSNB_REF.nc inputFiles/SnowMIP/cdp/cdp_init_run_${experiment}.nc
    cp inputFiles/SnowMIP/cdp/rsfile_run_TCZE_REF_TSNB_REF.nc inputFiles/SnowMIP/cdp/rsfile_run_${experiment}.nc

    # Launch simulation
    apptainer exec CLASSICcontainer.simg bin/CLASSIC_serial_TCZE_REF_TSNB_REF inputFiles/SnowMIP/cdp/job_options_file_run_${experiment}.txt 0/0 |& tee outputFiles/SnowMIP/cdp/run_${experiment}/output.txt
done

