# This script creates custom job options files for each fluxnet site in the inputFiles/FLUXNETsites directory

# Turn on Ncycle from command line.
# usage:
#  prep_jobopts.sh runfolder
#  prep_jobopts.sh runfolder -ncycle on

# Note: It is presently set to use the default shrub parameters file.

if [[ ($# -ne 1) && ($# -ne 3)]]; then
  echo ' usage for prep_jobopts.sh:'
  echo '    prep_jobopts.sh runfolder'
  echo '  or if you want to run with the N cycle on:'  
  echo '    prep_jobopts.sh runfolder -ncycle on'
  exit 2
fi 

Ncycle_toggle="off"

runfolder=$1

# Check if we are on the HPC. Use the uname, which will start with ppp if on the
# ECCC HPC. Otherwise, set false.
systest=`uname -n`
short=${systest:0:3}

if [[ "$short" == "ppp" ]]; then
  on_HPC=true
else
  on_HPC=false
fi

# Check if runfolder exists
if [[ ! -d $runfolder ]]; then
  echo 'stop: runfolder does not exist:'$runfolder
  exit 2
fi 

while :; do
  case $2 in
    -ncycle) 
             if [ "$3" ] ; then
               Ncycle_toggle=$3
               shift
             else
               Ncycle_toggle="off"
             fi
             ;;
    --)
        shift
        break
        ;;
    *)
        break
  esac
  shift
done

if [ "$Ncycle_toggle" = "on" ] ; then
  Ncycle=".true."
  deposition=".true."
else
  Ncycle=".false."
  deposition=".false."
fi
fertilizer=".false."
echo "Ncycle is set to: $Ncycle"

# First, we get a map of the repository by establishing where this script is located, then
# deducing where the root of the repository is.
script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
rootdir=${script_location%%/tools*}
inputs=$rootdir/inputFiles/FLUXNETsites
# We'll need the yaml reader to access the site_info.yaml for each site
yaml_reader=$script_location/yaml_reader.py
# Iterate through the site directories
for f in $inputs/*; do
  if [ -d $f ]; then
    # Get the name of the site (using its directory)
    current=${f##*/}
    echo "Setting up $current..."

    # Make the run folder for that site if not there already
    if [ ! -d $runfolder/$current ] ; then
      mkdir $runfolder/$current
    fi

    # Create the model init file from the cdl in the repo.
    if [[ ! -f $inputs/$current/$current\_init.cdl ]] ; then 
      echo 'Cannot find: '$inputs/$current/$current\_init.cdl 
    else 
      ncgen -o $runfolder/$current/$current\_init.nc  $inputs/$current/$current\_init.cdl 
      yes | cp $runfolder/$current/$current\_init.nc $runfolder/$current/rsfile.nc 
    fi 

    # Get metadata about the site for later use
    start_year=$(python3 $yaml_reader $f/siteinfo.yaml start)
    end_year=$(python3 $yaml_reader $f/siteinfo.yaml end)
    leap_flag=$(python3 $yaml_reader $f/siteinfo.yaml leap)
    ZRFH=$(python3 $yaml_reader $f/siteinfo.yaml ZRFH)
    ZRFM=$(python3 $yaml_reader $f/siteinfo.yaml ZRFM)

    jof=$runfolder/$current/job_options_file.txt
    echo "Copying template jobopts to $jof"
    yes | cp $rootdir/configurationFiles/template_job_options_file.txt $jof

    echo "Preparing $jof"
    # Now we use sed (stream editor) to tailor the job options file with the info we gained
    # from site_info.yaml
    # Replace start year and end year
    sed -i "/readMetStartYear/s| = [0-9]* | = $start_year | ; /readMetEndYear/s| = [0-9]* | = $end_year | ; /metLoop/s| = [0-9]* | = 1 |" $jof
    
    # Replace leap 
    sed -i "/leap/s| = \..*\, | = $leap_flag |" $jof

    # Use the proper default parameter file
    sed -i "/runparams_file/s|'.*'|'$rootdir/configurationFiles/default_run_parameters_shrubs.txt'|" $jof

    # Use the proper output variable descriptor file
    sed -i "/xmlFile/s|'.*'|'$rootdir/configurationFiles/outputVariableDescriptors.xml'|" $jof

    # Location of meteorological files is set based on whether or not the user is on HPC
    met_loc=$rootdir/inputFiles/$current
    if [ -n "$on_HPC" ] ; then
      met_loc="/space/hall5/sitestore/eccc/crd/ccrp/users/scrd530/benchmark_classic/CBC/sites/meteorology/${current}"
    fi
    
    # Replace all the meteorological files
    sed -i "/metFileFss/s|'.*'|'$met_loc/metVar_sw.nc'| ; /metFileFdl/s|'.*'|'$met_loc/metVar_lw.nc'| ; /metFilePre/s|'.*'|'$met_loc/metVar_pr.nc'| ; /metFileTa/s|'.*'|'$met_loc/metVar_ta.nc'| ; /metFileQa/s|'.*'|'$met_loc/metVar_qa.nc'| ; /metFileUv/s|'.*'|'$met_loc/metVar_wi.nc'| ; /metFilePres/s|'.*'|'$met_loc/metVar_ap.nc'|" $jof

    # Use the proper CO2 file
    if [ -n "$on_HPC" ] ; then
      sed -i "/CO2File/s|'.*'|'/space/hall5/sitestore/eccc/crd/ccrp/users/scrd530/classic_inputs/CO2/TRENDY_CO2_1700-2020_GCP2021.nc'|" $jof
    else
      sed -i "/CO2File/s|'.*'|'inputFiles/TRENDY_CO2_1700_2018.nc'|" $jof
    fi

    # Replace the init and restart files
    sed -i "/init_file/s|'.*'|'$runfolder/$current/${current}_init.nc'| ; /rs_file_to_overwrite/s|'.*'|'$runfolder/$current/rsfile.nc'|" $jof

    # Turn off fires and land use, and competition
    sed -i "/lnduseon/s|=\s*\..*\.\s*,|= \.false\. ,| ; /fixedYearLUC/s|=\s*[-0-9]*\s*,|= -9999 ,| ; /dofire/s|=\s*\..*\.\s*,|= \.false\. ,| ; /PFTCompetition/s|=\s*\..*\.|= \.false\.|" $jof

    # Replace the IDISP, IZREF, ZRFH, and ZRFM
    sed -i "/IDISP/s|= [0-9\.]*\s*,|= 1 ,| ; /IZREF/s|= [0-9\.]*\s*,|= 1 ,| ; /ZRFH/s|= [0-9\.]*,|= $ZRFH ,| ; /ZRFM/s|= [0-9\.]*,|= $ZRFM ,|" $jof
    
    # Replace Ncycle_on, fertilizeron, and depositionon
    sed -i "/Ncycle_on/s| = \..*\. | = $Ncycle | ; /fertilizeron/s| = \..*\. | = $fertilizer | ; /depositionon/s| = \..*\. | = $deposition |" $jof

    # Replace xmlFile
    sed -i "/xmlFile/s|'.*'|'$rootdir/configurationFiles/outputVariableDescriptors.xml'|" $jof
    
    # Replace output_directory
    sed -i "/output_director/s|'.*'|'$runfolder/outputFiles/$current/netCDF'|" $jof

    # Adjust Output Options to produce daily output only
    sed -i "/JDSTD/s|= [0-9\.]*\s*,|= 1 ,| ; /JDENDD/s|= [0-9\.]*\s*,|= 366 ,|" $jof
    sed -i "/JDSTY/s| = [0-9]* ,| = $start_year ,| ; /JDENDY/s| = [0-9]* ,| = $end_year ,|" $jof
    sed -i "/dodayoutput/s| = \..*\. | = .true. | ; /domonthoutput/s| = \..*\. | = .false. | ; /doAnnualOutput/s| = \..*\. | = .false. |" $jof

    echo "Making output directory"
    mkdir -p $runfolder/outputFiles/${current}/netCDF
  fi

done
