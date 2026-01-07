#!/bin/bash
set -e

# Script to launch CLASSIC model runs on ECCC collab (i.e. gpscc3).

# Usage: Modify the "User specified parameter" section as necessary and then just execute this script 
#        on the platform chosen to do the run.

# The final output netCDF files will be stored in the directory "netcdf_files" under the user specified
# output directory.

# -----------
# *** NB: AMBER/Edifice are not (yet?) implemented on collab.
# If a configure.R file exists in the user specified output directory, a job to run AMBER will be
# automatically launched on one of the ppps at the end of the model run.

# If an edifice.yaml file exists in the user specified output directory, a job to run Edifice will be
# automatically launched on one of the ppps at the end of the model run.
# -----------

# Background:

# CLASSIC currently runs well in parallel on a single node, but not across multiple nodes.

# The work-around is to run the model on a subset of the domain on each node. The number of nodes
# required is automatically computed based on the number of land grid cells in the domain, the
# estimated time required to run the simulation and the maximum wallclock allowed per node.

# At the end of the run, the files from each node will be stored in separate directories of the form 
# "netcdf_files_NN". Each directory will correspond to a subset of the full domain where "NN" is 
# numbered sequentially from 01, 02, etc. These are stitched together to produce the full domain 
# and will be stored in the directory "netcdf_files" under the user specified output directory.

# Initial implementation: Ed Chan, Feb 2024 (based on the script used for the ECCC HPC clusters)

# Modifications:

# ===============================================================================================
  
# -------------------------
# User specified parameters
# -------------------------

# Location of output files. 

root_dir=/space/hall0/work/eccc/crd/ccrp/$USER/classic-develop
output_directory=$root_dir/outputFiles/SnowArctic/tvc/run_1peat_2xSnowf_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne_correct_SH_collab

# Location of job options file. 
# *** NB: This script checks all input files (as specified in the job options file).
#         Files that are not needed for the run must either be commented out or specified as an empty string.

job_options_file=$root_dir/inputFiles/SnowArctic/tvc/job_options_file_run_1peat_2xSnowf_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne_correct_SH_collab.txt

# Location of compiled CLASSIC executable.
# *** NB: This script runs the model by first ensuring that the proper environment is set up.
#         The executable should be compiled using the same environment by running the script 
#         ~rec001/public/classic/make_classic.sh in the same directory as the Makefile.

executable=$root_dir/bin/CLASSIC_parallel_intel

# Specify the row/column limits OR the lat/lon limits of the domain/sub-domain.
# If not specified, defaults are set for the full domain.
# In the case of regular grids, the "full domain" is assumed to be global and
# ranges from latitudes -58 to 90 in order to avoid processing over Antarctica.

# Projected grids.
# This option is used only if the job options file contains the namelist parameter: projectedGrid=.true.
# Rows range from 1 to nlat and columns 1 to nlon.
# row1=45 ; row2=85 ; col1=157 ; col2=214

# Regular grids.
# Lons typically range from 0 to 360, but will depend on what is used in the input files. 
# Lats are typically ordered South to North, but will also depend the order used in the input files.
# lat1=-58 ; lat2=90 ; lon1=0 ; lon2=360 
# lat1=0 ; lat2=0 ; lon1=0 ; lon2=0 

# Email to send job start/end/abort notifications.

email="Mickael.Lalande@uqtr.ca"

# Job submission to gpscc3 allows the user to specify the size of the virtual file system "tmpfs" used
# to perform I/O while the model is running. Since tmpfs is allocated from the memory requested for each node 
# (i.e. 440 GB), the remainder is the amount of RAM available on each node. E.g. if tmpfs=200G, then the
# available RAM is: 440G - 200G = 240G. 
tmpfs=200G

# NB: Changing membal increases the number of nodes used for a fixed size of tmpfs. It may be better
#     to just increase the size of tmpfs as long as there is enough RAM remaining to run the model.
# Memory balance factor no more than (nrow/nodes)*membal will be assigned to any one node.
# If nodes are crashing due to memory issues membal can be set to between 2 and 3 to ensure balance
membal=2000

# ===============================================================================================

### ---------------------------------------------
### Lines below may require editing as necessary.
### ---------------------------------------------

# A metric "phi" is used below as an estimate of the model's performance on a given platform for a typical model configuration. 
# *** If jobs time out, it may need to be adjusted (downwards), especially if options such as saving daily output are turned on.
#
# The metric is defined as the number of grid cells that can be run over the length of the simulation 
# (in years) per hour (i.e. gridcells*years/hour) and can be estimated by timing a (short) model run. 
# It is used to estimate the number of nodes required to run the entire simulation under the constraint 
# of the max wallclock time permitted on each node.
#
# For the Canada domain, 18300 grid cells can run for 2 years in ~0.5 hours.
# If A = 18300, B = 2, and C = 0.5, then A*B/C = 73200 gridcell-years/hour.
# Choose a somewhat lower number to be safe.

#phi=65000
phi=60000

# Set parameters according to the platform.
# On gpscc3, wallclock is unlimited for the eccc_pegasus project, but the number of cores is limited to 768 (i.e. 12 full nodes). 
# However, for "low" priority (i.e. Quality of Service or QoS) jobs with wallclock times less than 6 hours,
# there is no restriction on the number of cores/nodes. This will be used as the default configuration.

[[ -z $JOBHOST ]] && JOBHOST=$ORDENV_TRUEHOST

case "$JOBHOST" in
  gpscc3) ncores=64
          max_wallclock=6 ;;
       *) echo "Current platform $JOBHOST is not supported. This script should only be used to submit jobs to gpscc3." && exit -1 ;;
esac

# -----------------------------------------------------------------------------------------------

### ------------------------------------------------
### Lines below should not normally require editing.
### ------------------------------------------------

# ---------------------------------------
# Check and/or prepare directories/files.
# ---------------------------------------

# Clear screen.

/usr/bin/clear

# Check if output netCDF directories exist and provide the option to either delete them or abort.

mkdir -p $output_directory
cd $output_directory

list=netcdf_files_[0-9]*
[ "$(echo $list)" = "$list" ] && unset list
[ -d netcdf_files ] && list="netcdf_files $list"

if [ -n "$list" ] ; then
  echo
  echo
  echo "The following directories in $output_directory will be deleted: "
  echo
  printf "%s\n" $list
  echo
  echo
  read -s -n1 -p "Press any key to continue or Ctrl-C to cancel ..."
  rm -rf $list
  echo
fi

# Set namelist parameters defined in the job options file as local shell variables. 
# NB: output_directory is defined in job options file and will override what is specified, so reset it.

new_out=$output_directory
eval "$(cat $job_options_file | sed -n -e 's/[&!]/#/g; s/,/;/g; s/ *= */=/gp')"
output_directory=$new_out

# Ensure all input files specified in the job options file are accessible. 
# *** NB: If a file is specified, it will be checked, even if it isn't actually used
#         in the model. In that case, a user can either specify an empty string or just
#         comment it out in the job options file. While it's possible to check 
#         switches to see if files are needed, this will be more involved as the 
#         current switch definitions don't contain common strings (like [Ff]file)
#         as for input files.

flag=false

for file_var in $(compgen -v | grep -G '[Ff]ile')
do
  infile=$(eval echo \$$file_var)
  [ -n "$infile" -a ! -f "$infile" ] && echo $infile is not accessible. && flag=true
done

$flag && echo 'Error: Job options file contains paths to inaccessible file(s)' && exit -2

# Ensure model executable is accessible. 

[ -n "$executable" -a ! -f "$executable" ] && echo $executable is not accessible. && exit -3

# Set the location of various files such as standard/error output from execution of the job, 
# copies of scripts used for the run, etc. Contents may be useful for debugging purposes.

run_files_directory=$output_directory/run_files
rm -rf $run_files_directory
mkdir $run_files_directory
cd $run_files_directory

# Set the job name to the name of the directory containing the run.

runname=${output_directory##*/}

# -----------------------------------------------------------------------------------------------

# --------------------------------------------------
# Obtain grid lat/lon and/or row/column information.
# --------------------------------------------------

# Extract the land mask from "init_file" to get grid-related information.
# Since cdo can't subset a "generic" grid, ensure that the grid type is set to "lonlat".

cdo -s -selvar,FLND $init_file GC.nc 2>/dev/null
cdo -s -griddes GC.nc | sed 's/generic/lonlat/' > griddes.txt
cdo -s -setgrid,griddes.txt GC.nc GC_tmp.nc
mv GC_tmp.nc GC.nc

# If the namelist parameter "projectedGrid" is set to ".true." in the job options file, 
# then following computations are based on rows/columns of a projected grid. Otherwise, 
# they are done based on lats/lons of a regular grid.

[ -n "$projectedGrid" ] && projectedGrid=$(echo $projectedGrid | tr '[:upper:]' '[:lower:]')

if [ "$projectedGrid" = '.true.' ] ; then

  [ -n "$lat1$lat2$lon1$lon2" ] && echo 'Subgrid limits lat1/lat2/lon1/lon2 not valid with projectedGrid = .true.' && exit -4

  # Get the number of rows/columns of the grid.

  nlat=$(ncdump -h GC.nc | fgrep "lat =" | cut -d' ' -f3)
  nlon=$(ncdump -h GC.nc | fgrep "lon =" | cut -d' ' -f3)

  # Set the row/column limits and/or defaults.

  [ -z "$row1" ] && row1=1 ; [ -z "$row2" ] && row2=$nlat ; [ $row2 -gt $nlat ] && row2=$nlat
  [ -z "$col1" ] && col1=1 ; [ -z "$col2" ] && col2=$nlon ; [ $col2 -gt $nlon ] && col2=$nlon

  # Calculate the number of land grid cells in the domain/sub-domain.

  total_land_cells=$(cdo -s -output -fldsum -setrtoc2,1.E-20,1.5,1,0 -setmisstoc,0 -selindexbox,$col1,$col2,$row1,$row2 GC.nc)

else

  [ -n "$row1$row2$col1$col2" ] && echo 'Subgrid limits row1/row2/col1/col2 not valid with projectedGrid = .false.' && exit -5

  # Set the lat/lon limits and/or defaults.

  lats=( $(ncks --cdl -C -v lat GC.nc | fgrep 'lat =' | tail -1 | sed -e 's/[;,]//g; s/^ *lat = //') )
  [ -z "$lat1" ] && lat1=${lats[0]} ; [ -z "$lat2" ] && lat2=${lats[-1]}
  [ $(echo "$lat1 < -58" | bc) == 1 ] && lat1=-58 && lat2=90
  [ $(echo "$lat2 < -58" | bc) == 1 ] && lat1=90 && lat2=-58

  lons=( $(ncks --cdl -C -v lon GC.nc | fgrep 'lon =' | tail -1 | sed -e 's/[;,]//g; s/^ *lon = //') )
  [ -z "$lon1" ] && lon1=${lons[0]}   ; [ -z "$lon2" ] && lon2=${lons[-1]}

  # Calculate the number of land grid cells in the domain/sub-domain.

  cdo -s -sellonlatbox,$lon1,$lon2,$lat1,$lat2 GC.nc GC_subdomain.nc
  total_land_cells=$(cdo -s -output -fldsum -setrtoc2,1.E-20,1.5,1,0 -setmisstoc,0 GC_subdomain.nc)

  # Get the row number (of the full domain) that corresponds to the first latitude of the sub-domain.

  lat1_subdomain=$(ncdump -v lat GC_subdomain.nc | fgrep 'lat =' | tail -1 | sed -e 's/,.*//; s/^ *lat = //')

  for i in ${!lats[@]}
  do
    [[ $lat1_subdomain == ${lats[$i]} ]] && row1=$((i+1))
  done

fi

# -----------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Determine the number of nodes, wallclock time, etc. required for the job.
# Also determine the number of grid rows to run on each node.
# -------------------------------------------------------------------------

# Compute the length of the run (years) based on information from the job options file.

run_length=$(echo "($readMetEndYear-$readMetStartYear+1)*$metLoop" | bc)

# Calculate the time required to do runs using the gridcell-years/hour metric.

total_time=$(echo "scale=2; $run_length*$total_land_cells/$phi" | bc)

# Estimate the number of nodes required. 

nodes=$(echo "$total_time/$max_wallclock" + 1 | bc)         # rounded up

# Loop to compute job parameters, giving the user an option to choose a different number of nodes.

flag=false

while :
do

  /usr/bin/clear
  
  # Get number of land grid cells in each row/lat of the grid. Note: cdo will only compute a zonal sum
  # for regular grids, so impose one on the data for projected grids, just for this computation.
  if [ "$projectedGrid" = '.true.' ] ; then
    land_cells_per_row=( $(cdo -s -output -zonsum -setrtoc2,1.E-20,1.5,1,0 -setmisstoc,0 -selindexbox,$col1,$col2,$row1,$row2 -setgrid,r${nlon}x$nlat GC.nc) )
  else
    land_cells_per_row=( $(cdo -s -output -zonsum -setrtoc2,1.E-20,1.5,1,0 -setmisstoc,0 -sellonlatbox,$lon1,$lon2,$lat1,$lat2 GC_subdomain.nc) )
  fi
  # make sure the number of nodes requested does not exceed the number of rows
  if [ $nodes -gt ${#land_cells_per_row[@]} ]; then nodes=$((${#land_cells_per_row[@]} - 1)); fi

  # Compute initial values of parameters used to estimate the amount of run time required on each node.

  grid_cells_per_node=$(echo "$total_land_cells/$nodes" | bc)
  run_time=$(echo "scale=2; $run_length*$grid_cells_per_node/$phi" | bc)
  
  # If the estimated run time is within 1/2 hour of the "max_wallclock", then it is safer to 
  # add an extra node and run fewer grid cells per node. 
  
  if (( $(echo "$max_wallclock-$run_time < 0.5" | bc -l) )) ; then
    grid_cells_per_node=$(echo "$total_land_cells/$((nodes+1))" | bc)
  fi
  
  # Calculate the maximum number of rows that can be assigned to a node without exceeding membal
  # the if ensures this value cannot be less than one
  rpnm=$(( (${#land_cells_per_row[@]}/$nodes)*$membal )); if [ $rpnm -lt 1 ]; then rpnm=1; fi

  # Get the row indices where the cumulative number of land grid cells between the indices do not exceed the
  # number of grid cells to be run per node. So, e.g. if the number of land grid cells in rows 1-50 are less 
  # than $grid_cells_per_node, then the 1st row index is 50. And if the same is true for rows 51-100, then
  # the next row index is 100. And if there are only the 2 indices, then rows 101-<last row> contains the 
  # last set of land grid cells (i.e. the job is split across 3 nodes).
  # the number of rows per node is also now not allowed to exceed rpnm to satisfy membal and prevent
  # overloaded nodes
  
  sum=0; row_boundaries=( $((row1-1)) ) ; land_cells_per_node=()  ; sum_row=0
  
  for i in ${!land_cells_per_row[@]}
  do
    sum0=$sum ; sum=$(( sum + land_cells_per_row[$i] )); sum_row=$(($sum_row + 1))
    if [ $sum -gt $grid_cells_per_node ] || [ $sum_row -gt $rpnm ]; then
      row_boundaries+=( $((i+row1-1)) )
      land_cells_per_node+=( $sum0 )
      sum=${land_cells_per_row[$i]}
      sum_row=1
    fi
  done
  
  row_boundaries+=($row2)
  land_cells_per_node+=($sum)
  nodes=${#land_cells_per_node[@]}

  # Now that the number of nodes is determined, compute the wallclock time per node for the job.
  # If the estimated run time is within 1/2 hour of the wallclock time, then it is safer to 
  # add an extra half-hour as long as the total doesn't exceed the "max_wallclock".
  
  grid_cells_per_node=$(echo "$total_land_cells/$nodes" | bc)
  half_hour_intervals=$(echo "$run_length*$grid_cells_per_node*2/$phi + 1" | bc) # rounded up
  
  if (( $(echo "$half_hour_intervals/2-$run_time < 0.5" | bc -l) && $half_hour_intervals < $max_wallclock*2 )) ; then
    half_hour_intervals=$((half_hour_intervals+1))
  fi
  
  run_time=$(echo "scale=2; $run_length*$grid_cells_per_node/$phi" | bc)
  wallclock=$( printf '%d:%02d:00\n' $(( half_hour_intervals/2 )) $(( (half_hour_intervals*30)%60 )) )

  # Output job parameters to the screen and obtain user feedback.
  
  echo "Total number of land grid cells: $total_land_cells"
  echo "Length of run (years): $run_length"
  echo
  echo "Total time required (hours): $total_time"
  echo "Number of nodes required: $nodes"
  echo "Estimated amount of time required per node (hours) = $run_time"
  echo "Amount of time requested per node (H:MM:SS) = $wallclock"
  echo "Maximum wallclock (hours): $max_wallclock"
  echo 
  read -p 'Choose one of the following options:

  - Enter/Return to continue with the above job configuration
  - Redo configuration by entering a different number of nodes to use (approximate) 
  - Ctrl-C to cancel

  Input: '

  [ $(echo "$run_time/$max_wallclock" | bc) -gt 0 ] && continue
  [ -z "$REPLY" ] && break

  while [[ ! $REPLY =~ ^[0-9]+$ || $REPLY = 0 ]]
  do 
    read -p "Only positive integers are allowed. Try again: "
  done

  nodes=$REPLY ; flag=true

done

# Allow user to also set the wallclock time requested per node to the max allowed, 
# if requesting a different number of nodes.

if $flag ; then
  echo
  read -p "Override wallclock time requested per node to the max of $max_wallclock hours (y or any other key to skip)?"
  if [[ $REPLY == y* ]] ; then
    wallclock="$max_wallclock:00:00"
    echo
    echo "Amount of time requested per node (H:MM:SS) = $wallclock"
  fi
fi

# -----------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------
# Construct the command lines used to run the model. Also save a table identifying the row/lat 
# boundaries associated with each node. These can be used to reconstruct the full domain.
# --------------------------------------------------------------------------------------------

# Since parameters passed to the model are different depending on the type of grid, the appropriate command line is constructed
# first and saved to a variable. Note that the -host option must be specified otherwise all commands will target
# the same node. The hosts can be obtained from $SLURM_NODELIST at runtime.

if [ "$projectedGrid" = '.true.' ] ; then

  printf "%-20s %-20s %-20s %-20s\n" 'Node Number' 'First Row' 'Last Row' 'Number of Grid Cells' > domain_partitions.txt

  for i in ${!land_cells_per_node[@]}
  do
    printf "%-20d %-20d %-20d %-20d\n" $((i+1)) $((row_boundaries[$i]+1)) ${row_boundaries[$((i+1))]} ${land_cells_per_node[$i]} >> domain_partitions.txt

    # Set the number of processes to the number of grid cells if less than $ncores or the SLURM job may hang.
    nprocs=$ncores
    [[ ${land_cells_per_node[$i]} -lt $ncores ]] && nprocs=${land_cells_per_node[$i]}

    command_lines+="time mpirun -n $nprocs -host \$host$i run_model_script $((i+1)) $col1 $col2 $((row_boundaries[$i]+1)) ${row_boundaries[$((i+1))]} &\n"
  done

else

  # Compute 1/2 the difference between each latitude (only valid for 2nd to last latitude).

  lats=( $lat1 ${lats[@]} )

  for i in ${!lats[@]}
  do
    delta_lat+=( $(echo "(${lats[$((i+1))]} - ${lats[$i]})/2" | bc -l) )
  done

  printf "%-20s %-20s %-20s %-20s\n" 'Node Number' 'First Lat' 'Last Lat' 'Number of Grid Cells' > domain_partitions.txt

  # Compute the boundaries for each latitude band and use these on the command lines. Using boundaries half-way between 
  # latitudes prevents the possibility of round-off errors causing problems.

  boundaries=( $(printf "%.4f" $lat1) )

  for i in ${!land_cells_per_node[@]}
  do

    if [ $i -eq $(( ${#land_cells_per_node[@]} - 1)) ] ; then
      boundaries+=( $(printf "%.4f" $lat2) )
    else
      boundaries+=( $(printf "%.4f" $(echo "${lats[${row_boundaries[$((i+1))]}]} + ${delta_lat[${row_boundaries[$((i+1))]}]}" | bc -l)) )
    fi

    printf "%-20d %-20.4f %-20.4f %-20d\n" $((i+1)) ${boundaries[$i]} ${boundaries[$((i+1))]} ${land_cells_per_node[$i]} >> domain_partitions.txt
    
    # Set the number of processes to the number of grid cells if less than $ncores or the SLURM job may hang.
    nprocs=$ncores
    [[ ${land_cells_per_node[$i]} -lt $ncores ]] && nprocs=${land_cells_per_node[$i]}

    command_lines+="time mpirun -n $nprocs -host \$host$i run_model_script $((i+1)) $lon1 $lon2 ${boundaries[$i]} ${boundaries[$((i+1))]} &\n"

  done

fi

command_lines=$(printf "$command_lines")

echo
cat domain_partitions.txt
echo

# -----------------------------------------------------------------------------------------------

# ----------------------------------------
# Construct the script and submit the job.
# ----------------------------------------

# Set copy utility for transferring files. 
# Use of mcp appears to give the best performance, but could use other utilities (e.g. sscp). 
# NB: mcp does not appear to be available on collab.
#copy=mcp
copy=cp

cat <<endjob > CLASSIC_${runname}.job
#!/bin/bash
#SBATCH --job-name=$runname
#SBATCH --partition=standard
#SBATCH --account=eccc_pegasus_crd
#SBATCH --qos=low
#SBATCH --output=$run_files_directory/${runname}_$(date +%F_%H-%M-%S)_${JOBHOST#*-}.out
#SBATCH --nodes=$nodes
#SBATCH --cpus-per-task=$ncores
#SBATCH --time=$wallclock
#SBATCH --mem=440G
#SBATCH --comment='tmpfs_size=$tmpfs'
#SBATCH --mail-user=$email
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

echo Job ID: \$SLURM_JOBID

# Load netCDF/HDF5 package, but first save a copy of LD_LIBRARY_PATH for later use with the stitcher.

. ssmuse-sh -x main/opt/intelcomp/inteloneapi-2023.2.0/intelcomp+mpi+mkl
libpath=\$LD_LIBRARY_PATH
. ssmuse-sh -x main/opt/hdf5-netcdf4/parallel/alllib/inteloneapi-2023.2.0/01

# Construct scripts to be run via MPI.

cd $run_files_directory

cat <<'endcat' > run_model_script
#!/bin/bash

# TMPDIR (tmpfs) on each compute node is automatically deleted when the job finishes. 
cd \$TMPDIR

if [ "\$MPI_LOCALRANKID" = 0 ] ; then

  # Only done for 1 MPI task per node.

  # Create run directory on TMPDIR and get local copy of restart file.
  # Also modify and save a copy of the job options file to use TMPDIR for the output_directory on each compute node.
  # NB: For each PBS job, TMPDIR is named the same on all compute nodes, but safer to evaluate TMPDIR during execution on each node.

  mkdir \$TMPDIR/$runname
  $copy $rs_file_to_overwrite $runname/rsFile_modified.nc
  sed "/output_directory/s|'.*'|'\$TMPDIR/$runname'| ; /rs_file_to_overwrite/s|'.*'|'\$TMPDIR/$runname/rsFile_modified.nc'|" $job_options_file > $runname/job_options.txt

  # NB: Work-around to mitigate slow I/O reads of met files: 
  #     1) Use ncks (with large chunk cache size) to extract only the lat/time slice required on each node into TMPDIR.
  #        Since tmpfs space is limited, using a Lempel-Ziv deflation/compression level of 1 greatly reduces the file sizes.
  #        Note that performance of CLASSIC on the compressed files can be impacted unless the chunk sizes (in lat/lon) are recduced to 1.
  #     2) Note that the -F option is used to enable Fortran indexing, which is required when the lats are rows (integers). 
  #        When the lats are reals (i.e. latitudes in degrees, the -F is ignored).
  #     3) Use of the --no_abc option is necessary, otherwise the order of the variables is changed and this causes problems with
  #        the canada domain!?
  #     4) Modify the job options file accordingly.

  ifiles='$metFileFss $metFilefracFsf $metFileFdl $metFilePre $metFileSnow $metFileTa $metFileQa $metFileUv $metFilePres' 
  for ifile in \$ifiles ; do
    ofile=\${ifile##*/}
    ncks -h -F --no_abc -L1 --cnk_csh=1000000000 --cnk_map=dmn --cnk_plc=g3d --cnk_dmn=lat,1 --cnk_dmn=lon,1 \
         -d lat,\$4,\$5 -d time,${readMetStartYear}0000.0,${readMetEndYear}1232.0 \$ifile \$ofile
    sed "s|\$ifile|\$TMPDIR/\$ofile|" $runname/job_options.txt > $runname/opts.txt
    mv $runname/opts.txt $runname/job_options.txt
  done

  # Get listing of met files.
  #ls -alh

fi

# Run the model (all MPI tasks per node). 

$executable $runname/job_options.txt \$2/\$3/\$4/\$5

# Copy files to the temporary run directory (1 MPI task per node).
# No need to delete files on tmpfs (should be done by the system).

if [ "\$MPI_LOCALRANKID" = 0 ] ; then
  run_directory=netcdf_files_\$(echo 0\$1 | tail -c3)
  mv $runname \$run_directory
  $copy -r \$run_directory $output_directory
fi
endcat

chmod u+x run_model_script

# Distribute the execution of the model across all nodes. 

node_range=\${SLURM_NODELIST#*-}
node1=\$(echo \$node_range | sed 's/\[//; s/-.*//')
node1=\$((10#\$node1)) # remove leading zeros otherwise variable is interpreted as a string
cluster=\${SLURM_NODELIST%%-*}

for ((i = 0 ; i < \$SLURM_NNODES ; i++)); do
  node=\$(( node1 + i ))
  node=\$(echo 00\$node | tail -c4)
  eval host\$i=\$cluster-\$node
  eval echo \\\$host\$i
done

echo 'Executing CLASSIC'
$command_lines
wait

# Save a copy of the modified job_options file.

cp $output_directory/netcdf_files_01/job_options.txt $run_files_directory

# If necessary, stitch row/lat bands split across multiple output directories onto the full domain.
# NB: cdo/nco utilities used by the stitcher are not compatible with the parallel netCDF/HDF libaries 
#     used to compile/run the model. So reset the path first.

echo 'Executing classic_stitch_netcdf.sh'
export LD_LIBRARY_PATH=\$libpath
time ~rec001/public/classic/classic_stitch_netcdf.sh $output_directory 32

# If a configure.R file is present in the output directory, then submit job to run AMBER.

#front_end=\$ORDENV_TRUEHOST
#[ \$ORDENV_TRUEHOST == underhill ] && front_end=ppp5
#[ \$ORDENV_TRUEHOST == robert ] && front_end=ppp6

#if [ -f $output_directory/configure.R ] ; then
#  ~rec001/public/classic/amber.sh -c \$front_end $output_directory
#fi

# If an edifice.yaml file is present in the output directory, then submit job to run Edifice.
#if [ -f $output_directory/edifice.yaml ] ; then
#  ~rec001/public/classic/edifice.sh -c \$front_end $output_directory
#fi

endjob

# Clean up.

rm -f GC.nc GC_subdomain.nc griddes.txt

# Submit the job.

jobsub -c $JOBHOST CLASSIC_${runname}.job
   
echo
echo 'Job submission is complete.'
echo 
echo 'Check job status with:' 
echo '~rec001/public/classic/jobstat'
echo
echo 'Check job output files for model text and error messages in:'
echo $run_files_directory
echo

exit
