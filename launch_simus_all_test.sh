#!/bin/bash

# List of experiments
experiments=(
    #PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne # best with TVC
    #PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne_test_4
    #ALL_TVC_VA_MOST
    ALL_TVC_fact
    ALL_TVC_ZRF_fix
   
)

# List of sites
sites=(
    #"SnowMIP/cdp"
    #"SnowMIP/rme"
    #"SnowMIP/snb"
    #"SnowMIP/swa"
    #"SnowMIP/sap"
    #"SnowMIP/sod"
    #"SnowMIP/wfj"
    #"SnowArctic/byl"
    "SnowArctic/umt"
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
        if [ "${site}" == "SnowArctic/umt" ]; then
            cp inputFiles/${site}/job_options_file_${run_prefix}_Ref_30min_ext_correct.txt inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt
        else
            cp inputFiles/${site}/job_options_file_${run_prefix}_Ref_30min_ext.txt inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt
        fi

        # Modify path in job_options_file
        sed -i "s/Ref_30min_ext/${experiment}/g" inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt
        sed -i "s/dohhoutput = .true./dohhoutput = .false./g" inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt
        sed -i "s/domonthoutput = .true./domonthoutput = .false./g" inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt
        sed -i "s/doAnnualOutput = .true./doAnnualOutput = .false./g" inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt

        # Create output file directory
        mkdir -p outputFiles/${site}/${run_prefix}_${experiment}

        # Copy and restart files from REF run
        cp inputFiles/${site}/${site_base}_init_${run_prefix}_Ref_30min_ext.nc inputFiles/${site}/${site_base}_init_${run_prefix}_${experiment}.nc
        cp inputFiles/${site}/rsfile_${run_prefix}_Ref_30min_ext.nc inputFiles/${site}/rsfile_${run_prefix}_${experiment}.nc

        # Launch simulation
        apptainer exec CLASSICcontainer.simg bin/CLASSIC_serial_${experiment} inputFiles/${site}/job_options_file_${run_prefix}_${experiment}.txt 0/0 |& tee outputFiles/${site}/${run_prefix}_${experiment}/output.txt
    done
done


