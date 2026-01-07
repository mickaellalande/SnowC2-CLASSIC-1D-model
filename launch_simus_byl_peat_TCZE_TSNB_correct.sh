#!/bin/bash

# List of experiments
list=(
    TCZE_OP1_TSNB_OP5
)

# Loop over each experiment in the list
for experiment in "${list[@]}"; do
    # Print the current experiment
    echo "Processing experiment: ${experiment}"
    
    # Copy the job_option_file from ref run
    cp inputFiles/SnowArctic/byl/job_options_file_run_peat_TCZE_REF_TSNB_REF.txt inputFiles/SnowArctic/byl/job_options_file_run_peat_${experiment}.txt

    # Modify path in job_options_file
    sed -i "s/TCZE_REF_TSNB_REF/${experiment}/g" inputFiles/SnowArctic/byl/job_options_file_run_peat_${experiment}.txt

    # Create output file directory
    mkdir -p outputFiles/SnowArctic/byl/run_peat_${experiment}

    # Copy and restart files from REF run
    cp inputFiles/SnowArctic/byl/byl_init_run_peat_TCZE_REF_TSNB_REF.nc inputFiles/SnowArctic/byl/byl_init_run_peat_${experiment}.nc
    cp inputFiles/SnowArctic/byl/rsfile_run_peat_TCZE_REF_TSNB_REF.nc inputFiles/SnowArctic/byl/rsfile_run_peat_${experiment}.nc

    # Launch simulation
    apptainer exec CLASSICcontainer.simg bin/CLASSIC_serial_${experiment} inputFiles/SnowArctic/byl/job_options_file_run_peat_${experiment}.txt 0/0 |& tee outputFiles/SnowArctic/byl/run_peat_${experiment}/output.txt
done

