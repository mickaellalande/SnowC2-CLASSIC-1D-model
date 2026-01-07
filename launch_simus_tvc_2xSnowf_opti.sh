#!/bin/bash

# List of experiments
experiments=(
    #PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_3.0_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_5.0_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_4.0_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_3.5_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_4.0_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.0_calonne
    #DEF
    #BUG_CORRECT
    #BUG_CORRECT_TSNBT_OP1
    #BUG_CORRECT_TSNBT_OP1_EZERO
    #PHYS_ALL_SUBLI_v2
    #PHYS_ALL_SUBLI_v2_COMPAC_v1
    PHYS_ALL_SUBLI_CORRECT
    PHYS_ALL_SUBLI_CORRECT_COMPAC
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
        cp inputFiles/${site}/job_options_file_${run_prefix}_2xSnowf.txt inputFiles/${site}/job_options_file_${run_prefix}_2xSnowf_${experiment}.txt

        # Modify path in job_options_file
        sed -i "s/1peat_2xSnowf/1peat_2xSnowf_${experiment}/g" inputFiles/${site}/job_options_file_${run_prefix}_2xSnowf_${experiment}.txt
		
        # Create output file directory
        mkdir -p outputFiles/${site}/${run_prefix}_2xSnowf_${experiment}

        # Copy and restart files from REF run
        cp inputFiles/${site}/${site_base}_init_${run_prefix}_30min_ext.nc inputFiles/${site}/${site_base}_init_${run_prefix}_2xSnowf_${experiment}.nc
        cp inputFiles/${site}/rsfile_${run_prefix}_30min_ext.nc inputFiles/${site}/rsfile_${run_prefix}_2xSnowf_${experiment}.nc

        # Launch simulation
        apptainer exec CLASSICcontainer.simg bin/CLASSIC_serial_${experiment} inputFiles/${site}/job_options_file_${run_prefix}_2xSnowf_${experiment}.txt 0/0 |& tee outputFiles/${site}/${run_prefix}_2xSnowf_${experiment}/output.txt
    done
done


