# This runs over all the sites using CLASSIC. It can either use the container to run them in, or just use the 
# CLASSIC executable alone. 
# Pass the same run folder given to prep_jobopts.sh and the folder the container will use as its root folder or the executable path.
# If using a container, note that the root folder for the container must have both your CLASSIC folder and run folder in it so the container can see them both.
# Usage:
#   run_sites.sh runfolder container_root/executable path

runfolder=$1
container_root=$2

# Set default
useContainer=true

# Set the executable
exec=CLASSIC_serial_intel
#exec=CLASSIC_serial_gnu

# Check if we are running with the container (as in one was given)
if [[ ($# -ne 1) && ($# -ne 2)]]; then
  echo ' Usage: run_sites.sh runfolder container_root/executable path'
  echo ' or, exclude the container path if not wanting to use a container'
  exit 2
fi 
if [[ ($# -eq 1)]]; then
  echo 'No container given, so preparing run without one.'
  useContainer=false
fi 

# Check if runfolder exists
if [[ ! -d $runfolder ]]; then
  echo 'stop: runfolder does not exist:'$runfolder
  exit 2
fi

if [[ (! -d $container_root) && ("$useContainer" == "true") ]]; then
  echo 'stop: the folder you would like to use as the root does not exist:'$container_root
  exit 2
fi

script_location="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
rootdir=${script_location%%/tools*}
container=$rootdir/CLASSIC_container.simg

# Iterate through the sites and start a run for each one sequentially (using the container)
for f in $runfolder/*; do
  current=${f##*/}
  if [[ -d $f && ${current} != "outputFiles" ]]; then
    echo
    echo
    echo "Running $current..."
    echo
    if [[ "$useContainer" == "true" ]]; then
      singularity exec -B $container_root:$container_root $container $rootdir/bin/$exec $f/job_options_file.txt 0/0
    else
      $rootdir/bin/$exec $f/job_options_file.txt 0/0
    fi
  fi
done
