#!/bin/bash

# List of experiments
experiments=(
    DEF
    BUG_CORRECT
    BUG_CORRECT_TSNBT_OP1
    BUG_CORRECT_TSNBT_OP1_EZERO
    PHYS_ALL_SUBLI_v2
    PHYS_ALL_SUBLI_v2_COMPAC_v1
    PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne # COMPACT_2.5 best (without TVC)
    PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne # best with TVC
)

# List of sites
sites=(
    "SnowArctic/umt"
)

# Loop over each site
for site in "${sites[@]}"; do
    # Loop over each experiment in the list
    for experiment in "${experiments[@]}"; do
        # Define modified experiment name
        exp_mod="${experiment}_correct_SH"

        # Print the current site and experiment
        echo "Processing site: ${site}, experiment: ${exp_mod}"

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
            cp inputFiles/${site}/job_options_file_${run_prefix}_Ref_30min_ext_correct.txt inputFiles/${site}/job_options_file_${run_prefix}_${exp_mod}.txt
        else
            cp inputFiles/${site}/job_options_file_${run_prefix}_Ref_30min_ext.txt inputFiles/${site}/job_options_file_${run_prefix}_${exp_mod}.txt
        fi

        # Modify path in job_options_file
        sed -i "s/Ref_30min_ext/${exp_mod}/g" inputFiles/${site}/job_options_file_${run_prefix}_${exp_mod}.txt
        sed -i "s/dohhoutput = .true./dohhoutput = .false./g" inputFiles/${site}/job_options_file_${run_prefix}_${exp_mod}.txt
        sed -i "s/domonthoutput = .true./domonthoutput = .false./g" inputFiles/${site}/job_options_file_${run_prefix}_${exp_mod}.txt
        sed -i "s/doAnnualOutput = .true./doAnnualOutput = .false./g" inputFiles/${site}/job_options_file_${run_prefix}_${exp_mod}.txt

        # Modify forcing in job_options_file
        sed -i "s/met_insitu_umt_30min_ext_2012_2021_Qair/met_insitu_umt_30min_ext_2012_2021_Qair_correct_SH/g" inputFiles/${site}/job_options_file_${run_prefix}_${exp_mod}.txt

        # Create output file directory
        mkdir -p outputFiles/${site}/${run_prefix}_${exp_mod}

        # Copy and restart files from REF run
        cp inputFiles/${site}/${site_base}_init_${run_prefix}_Ref_30min_ext.nc inputFiles/${site}/${site_base}_init_${run_prefix}_${exp_mod}.nc
        cp inputFiles/${site}/rsfile_${run_prefix}_Ref_30min_ext.nc inputFiles/${site}/rsfile_${run_prefix}_${exp_mod}.nc

        # Launch simulation
        apptainer exec CLASSICcontainer.simg bin/CLASSIC_serial_${experiment} inputFiles/${site}/job_options_file_${run_prefix}_${exp_mod}.txt 0/0 |& tee outputFiles/${site}/${run_prefix}_${exp_mod}/output.txt
    done
done

