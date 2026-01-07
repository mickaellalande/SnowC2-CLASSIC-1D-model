#!/bin/bash

# List of experiments
experiments=(
    #RHOMAX_350
    #RHOMAX_700
    #RHOMAX_WIND_10
    #RHOMAX_WIND_5
    #RHO_WIND_THO_2
    #RHOMAX_WIND_5_THO_2
    #RHOSNI_CROCUS_RHOMAX_WIND_5
    #RHOSNI_CROCUS_R21_RHOMAX_WIND_5
    #RHOSNI_CROCUS_RHOMAX_WIND_5_TSNB_OP1
    #RHOSNI_CROCUS_RHOMAX_WIND_5_TSNB_OP1_EZERO
    #RHOSNI_CROCUS_RHOMAX_WIND_5_TSNB_OP1_EZERO_SUBLI
    #RHOSNI_CROCUS_RHOMAX_WIND_4_TSNB_OP1_EZERO_SUBLI
    #RHOSNI_CROCUS_RHOMAX_WIND_3_TSNB_OP1_EZERO_SUBLI
    #RHOSNI_CROCUS_RHOMAX_WIND_2_TSNB_OP1_EZERO_SUBLI
    #RHOSNI_CROCUS_RHOMAX_WIND_1_TSNB_OP1_EZERO_SUBLI
    #RHOSNI_CROCUS_RHOMAX_WIND_2.5_TSNB_OP1_EZERO_SUBLI
    #RHOMAX_350_ALL
    #RHOMAX_350_SD_ALL
    #RHOMAX_450_SD_ALL
    #RHOMAX_GAUSS_0.5_SIGMA_0.5_WIND_5_ALL
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL
    #RHOMAX_GAUSS_0.2_SIGMA_0.5_WIND_2.5_ALL
    #RHOMAX_GAUSS_0.2_SIGMA_0.2_WIND_2.5_ALL
    #RHOSNI_CROCUS_RHOMAX_WIND_2.5_TSNB_OP1_EZERO_1_SUBLI
    #RHOSNI_CROCUS_RHOMAX_WIND_2.5_TSNB_OP1_EZERO_SUBLI_Calonne
    #RHOSNI_CROCUS_RHOMAX_WIND_2.5_TSNB_OP1_EZERO_SUBLI_nofreeze
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_9_RATE_3
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_9_RATE_3_COMPAC_600
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_5_RATE_3_COMPAC_600
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_5_RATE_3_COMPAC_450
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_5_RATE_2_COMPAC_450
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.4_ALL_LIMIT_2.4
    RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_400
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


