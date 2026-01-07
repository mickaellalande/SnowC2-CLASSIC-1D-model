#!/bin/bash

# List of experiments
list=(
    TCZE_REF_TSNB_OP1
    TCZE_REF_TSNB_OP2
    TCZE_REF_TSNB_OP3
    TCZE_REF_TSNB_OP4
    TCZE_REF_TSNB_OP5
    TCZE_OP1_TSNB_REF
    TCZE_OP1_TSNB_OP1
    TCZE_OP1_TSNB_OP2
    TCZE_OP1_TSNB_OP3
    TCZE_OP1_TSNB_OP4
    TCZE_OP1_TSNB_OP5
    TCZE_OP2_TSNB_REF
    TCZE_OP2_TSNB_OP1
    TCZE_OP2_TSNB_OP2
    TCZE_OP2_TSNB_OP3
    TCZE_OP2_TSNB_OP4
    TCZE_OP2_TSNB_OP5
    TCZE_OP3_TSNB_REF
    TCZE_OP3_TSNB_OP1
    TCZE_OP3_TSNB_OP2
    TCZE_OP3_TSNB_OP3
    TCZE_OP3_TSNB_OP4
    TCZE_OP3_TSNB_OP5
    TCZE_REF_TSNB_REF_calonne2011
)

# Loop over each experiment in the list
for experiment in "${list[@]}"; do
    # Print the current experiment
    echo "Processing experiment: ${experiment}"
    
    # Copy the job_option_file from ref run
    cp inputFiles/SnowMIP/cdp/job_options_file_run_TCZE_REF_TSNB_REF.txt inputFiles/SnowMIP/cdp/job_options_file_run_${experiment}.txt

    # Modify path in job_options_file
    sed -i "s/TCZE_REF_TSNB_REF/${experiment}/g" inputFiles/SnowMIP/cdp/job_options_file_run_${experiment}.txt

    # Create output file directory
    mkdir -p outputFiles/SnowMIP/cdp/run_${experiment}

    # Copy and restart files from REF run
    cp inputFiles/SnowMIP/cdp/cdp_init_run_TCZE_REF_TSNB_REF.nc inputFiles/SnowMIP/cdp/cdp_init_run_${experiment}.nc
    cp inputFiles/SnowMIP/cdp/rsfile_run_TCZE_REF_TSNB_REF.nc inputFiles/SnowMIP/cdp/rsfile_run_${experiment}.nc

    # Launch simulation
    apptainer exec CLASSICcontainer.simg bin/CLASSIC_serial_${experiment} inputFiles/SnowMIP/cdp/job_options_file_run_${experiment}.txt 0/0 |& tee outputFiles/SnowMIP/cdp/run_${experiment}/output.txt
done

