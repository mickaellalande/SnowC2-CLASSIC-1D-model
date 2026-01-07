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
    #RHOMAX_GAUSS_0.5_WIND_2.5_ALL
    #RHOMAX_GAUSS_0.5_WIND_2.5_MAXWET_600_ALL
    #RHOMAX_GAUSS_0.5_SIGMA_0.5_WIND_5_ALL
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_3_ALL
    #RHOMAX_GAUSS_0.2_SIGMA_0.5_WIND_2.5_ALL
    #RHOMAX_GAUSS_0.2_SIGMA_0.2_WIND_2.5_ALL
    #RHOSNI_CROCUS_RHOMAX_WIND_2.5_TSNB_OP1_EZERO_1_SUBLI
    #RHOMAX_GAUSS_0.5_WIND_2.5_EZERO_1_ALL
    #RHOSNI_CROCUS
    #RHOSNI_CROCUS_TSNB_OP1
    #RHOSNI_CROCUS_TSNB_OP1_EZERO
    #RHOSNI_CROCUS_TSNB_OP1_EZERO_SUBLI
    #RHOSNI_W24
    #RHOMAX_GAUSS_0.5_WIND_2.5_ALL_LIMIT_3
    #RHOMAX_GAUSS_0.5_WIND_2.5_ALL_CORRECT
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_3_ALL_LIMIT_3
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_3_ALL_LIMIT_3_W24
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_5_RATE_2_COMPAC_450
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.6_ALL_LIMIT_2.6
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.7_ALL_LIMIT_2.7
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.8_ALL_LIMIT_2.8
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_400_calonne
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_calonne
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_400
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_425_calonne
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_calonne_R21
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_425_calonne_R21
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_425_calonne_GW1
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_425_calonne_GW2
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_425_calonne_W24
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_425_calonne_DEF
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_410_calonne_DEF
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_420_calonne_DEF
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_0.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_440_calonne_DEF
    #RHOMAX_GAUSS_0.5_SIGMA_0.8_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_0.5_SIGMA_0.9_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_0.5_SIGMA_1.1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_0.5_SIGMA_1.2_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_0.3_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_0.4_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_0.6_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_0.7_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_0.8_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_0.9_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_1.0_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_1.5_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #RHOMAX_GAUSS_2.0_SIGMA_1_WIND_2.5_ALL_LIMIT_2.5_COMPAC_430_calonne_DEF
    #DEF
    #BUG_CORRECT
    #BUG_CORRECT_TSNBT_OP1
    #BUG_CORRECT_TSNBT_OP1_EZERO
    #PHYS_ALL_SUBLI
    #PHYS_ALL_SUBLI_COMPAC
    #PHYS_ALL_SUBLI_COMPAC_calonne
    #PHYS_ALL_calonne
    #PHYS_ALL_TCZERO
    #PHYS_ALL_TCZERO_calonne
    #PHYS_ALL_TCZERO_calonne_600
    #PHYS_ALL_SUBLI_CORRECT_COMPAC_calonne
    #SUBLI_test
    #RHOSNI_SnowTran
    #RHOSNI_SnowTran_Arctic
    #SCF_3
    #PHYS_ALL_SUBLI_CORRECT
    #PHYS_ALL_SUBLI_CORRECT_COMPAC
    #PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne # COMPACT_2.5 best (without TVC)
    #PHYS_ALL_SUBLI_v2_COMPAC_4.0_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_3.0_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_3.5_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne # best with TVC
    #PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_4.0_calonne
    #PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.0_calonne
    PHYS_ALL_SUBLI_v2
    PHYS_ALL_SUBLI_v2_COMPAC_v1
   
)

# List of sites
sites=(
    "SnowMIP/cdp"
    "SnowMIP/rme"
    "SnowMIP/snb"
    "SnowMIP/swa"
    "SnowMIP/sap"
    "SnowMIP/sod"
    "SnowMIP/wfj"
    "SnowArctic/byl"
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


