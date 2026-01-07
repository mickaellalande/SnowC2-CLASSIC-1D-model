#!/bin/bash
# Pass the same run folder given to prep_jobopts.sh and run_sites.sh as well as the folder the container will use as its root folder.
# Note that the root folder for the container must have both your CLASSIC folder and run folder in it so the container can see them both.
# Usage:
#   process_outputs.sh runfolder(s)container_root (if using one)
# Note: the plots are placed into the first specified runfolder

# Find the run folders by making an array of the arguments
runfolder=( "$@" )
#printf "%s\n" "${runfolder[@]}"

# Container would be the last one.
container_root=${@: -1}
a=$(echo $container_root | awk '/.simg/{i++}i==1{print; exit}')
if [ -z $a ]; then
  useContainer=false
else
  useContainer=true
  # remove the container name from the list of runfolders
  unset runfolder[-1] 
fi

# List all to plot here. Should also be in the plot_site notebook.
vars_to_plot=('Latent_Heat' 'Sensible_Heat' 'Gross_Primary_Productivity' 'Ecosystem_Respiration' 'Net_Ecosystem_Productivity' 'Leaf_Area_Index' 'Precipitation' 'Air_Temperature')

# Specify where the obs are located. It tries them in order and stops at the first place it finds the needed site file.
obs_path="'/space/hall5/sitestore/eccc/crd/ccrp/users/scrd530/benchmark_classic/CBC/sites/daily_observations/FLUXNET2015_daily_obs_all_sites', '/space/hall5/sitestore/eccc/crd/ccrp/users/scrd530/benchmark_classic/CBC/sites/daily_observations/AmeriFlux_beta_May2022_daily_obs', '/space/hall5/sitestore/eccc/crd/ccrp/users/scrd530/benchmark_classic/CBC/sites/daily_observations/Collaborator_sites_daily_obs'"

# Path to conda init
ssmuse=". ssmuse-sh -x /fs/ssm/eccc/crd/ccrp/ssm/miniconda3_22.9.0_all"

# Conda environment to use
condaenv="edifice_core"

if [[ ($# -eq 0) ]]; then
  echo ' Usage: process_outputs.sh runfolder(s) container_root'
  echo ' or, exclude the container path if not wanting to use a container'
  echo ' container_root is always specified last'
  exit 2
fi 

# Check if runfolder(s) exist(s)
for i in "${runfolder[@]}" ; do
  if [[ ! -d $i ]]; then
    echo 'stop: runfolder does not exist:'$runfolder
    exit 2
  fi
done

if [[ ! -d $container_root && "$useContainer" == "true" ]]; then
  echo 'stop: the folder you would like to use as the root for the container does not exist:'$container_root
  exit 2
fi

script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
rootdir=${script_location%%/tools*}
container=/home/jskye99/ECCCmaterial/CLASSIC-developFLUXNETupdate/tools/singularityContainerRecipe/newpackagetest.simg #$rootdir/CLASSIC_container.simg

if [[ "$useContainer" == "false" ]]; then # find the conda and activate the edifice environment
  # Enter the conda environment of your choice:
  $ssmuse
  conda activate $condaenv
fi

# Iterate through the FLUXNET output directories. Any non-empty directories will
# have their outputs converted to csv files, then moved into an appropriate csv
# sub-directory
for i in "${runfolder[@]}" ; do
  for f in $i/outputFiles/*; do
    rm -rf $f/csv
    files=$(ls $f/netCDF | wc -l)
    [ -z "$(ls $f/netCDF | grep .nc)" ] && continue
    #echo "Converting netCDF for ${f##*/}"
    echo "Converting netCDF for $f"
    if [ "$files" -gt "0" ]; then
      if [[ "$useContainer" == "true" ]]; then
        singularity exec -B $container_root:$container_root $container python3 $rootdir/tools/convertOutputToCSV/convertNetCDFOutputToCSV_batch.py $f/netCDF
      else
        python3 $rootdir/tools/convertOutputToCSV/convertNetCDFOutputToCSV_batch.py $f/netCDF
      fi
      mkdir -p $f/csv
      yes | mv $f/netCDF/*.csv $f/csv
    fi
  done
done


# # Run the plot generator on the FLUXNET outputs. Output may be quite substantial
# # if not all sites have output. This is to be expected and is not a problem.
# echo "Creating folder for plots..."
# We place the plots_dir in the folder of the first specified run by default.

plots_dir=${runfolder[0]}/outputFiles/plots

# # plots_dir=/home/rjs001/public_html/quick_start_test/plots
rm -rf $plots_dir
mkdir $plots_dir

echo "Creating world maps of site locations..."
if [[ "$useContainer" == "true" ]]; then
  singularity exec -B $container_root:$container_root $container python3 $rootdir/tools/siteLevelFLUXNET/worldmap.py $rootdir/inputFiles/FLUXNETsites $runfolder/outputFiles/plots
else
  python3 $rootdir/tools/siteLevelFLUXNET/worldmap.py $rootdir/inputFiles/FLUXNETsites $runfolder/outputFiles/plots
fi

echo "Creating Canada maps of site locations..."
if [[ "$useContainer" == "true" ]]; then
  singularity exec -B $container_root:$container_root $container python3 $rootdir/tools/siteLevelFLUXNET/canadamap.py $rootdir/inputFiles/FLUXNETsites $runfolder/outputFiles/plots
else
  python3 $rootdir/tools/siteLevelFLUXNET/canadamap.py $rootdir/inputFiles/FLUXNETsites $runfolder/outputFiles/plots
fi 

echo "Creating Timeseries, Crossplots, and Seasonal Coverage graphs..."
ntb=$rootdir/tools/siteLevelFLUXNET/plot_site.ipynb

#sed -i "/obs_data_path =/s|'.*'|'$rootdir/inputFiles/obsFiles'|" $ntb
sed -i "/obs_data_path =/s|'.*'|$obs_path|" $ntb
sed -i "/yaml_data_path =/s|'.*'|'$rootdir/inputFiles/FLUXNETsites'|" $ntb

# Use the first runfolder as the entry point
for f in ${runfolder[0]}/outputFiles/*; do
  current=${f##*/}
  if [ ${current} == 'plots' ]; then
    continue
  fi
  rm -rf $plots_dir/$current
  mkdir $plots_dir/$current
  ntb_current=$plots_dir/$current/plot_site_$current.ipynb
  cp $ntb $ntb_current
  # Now add in the data as lists for the different runs, if more than one.
  if (( ${#runfolder[@]} > 1 )); then
    for i in "${runfolder[@]}" ; do
      run=$(basename $i)
      namz="$namz, '$run'"
      filez="$filez, '$i/outputFiles/$current/csv'"
    done 
    # remove the first comma
    namz="${namz:1}"
    filez="${filez:1}"
  else
    filez="'"$f/csv"'"
    namz="'"$(basename ${runfolder[0]})"'"
  fi 
  sed -i "/model_data_path =/s|'.*'|$filez|" $ntb_current
  sed -i "/runnames =/s|'.*'|$namz|" $ntb_current
  #sed -i "/model_data_path =/s|'.*'|'$f/csv'|" $ntb_current
  #sed -i "/runnames =/s|'.*'|'$f/csv'|" $ntb_current
  sed -i "/site_name =/s|'.*'|'$current'|" $ntb_current
  # exit
  for v in "${vars_to_plot[@]}"; do
    sed -i "/var_to_plot =/s|'.*'|'$v'|" $ntb_current
    if [[ "$useContainer" == "true" ]]; then
      singularity exec -B $container_root:$container_root $container jupyter nbconvert --execute --no-input --to html --output=${current}_$v.html $ntb_current
    else 
      jupyter nbconvert --execute --no-input --to html --output=${current}_$v.html $ntb_current
    fi
  done
done

cd $plots_dir
if [[ "$useContainer" == "true" ]]; then
  singularity exec -B $container_root:$container_root $container tree -H "." -L 1 --noreport --charset utf-8 -P "*.html" -o $plots_dir/index.html
else 
  tree -H "." -L 1 --noreport --charset utf-8 -o $plots_dir/index.html
fi
