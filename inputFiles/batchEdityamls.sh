# This script runs through all of the FLUXNET2015 sites we have in the CBC
# and edits their yaml files

# J Melton May 2022

#cbcsites=/home/rjm001/CLASSIC/inputFiles/FLUXNETsites
cbcsites=/home/rjm001/CLASSIC/inputFiles/FLUXNETsites_in_progress/COHERENT_Cinit

# Iterate through the site directories
for f in $cbcsites/*; do
  if [ -d $f ]; then
    # Get the name of the site (using its directory)
    current=${f##*/}
    echo $current 
    # put 3 new lines after leap for peatland, landscape_state, doi
    sed -i '/^leap.*/a peatland: "false"\nlandscape_state: "undisturbed"\ndoi: "10.3334/ORNLDAAC/1335"' $f/siteinfo.yaml 

  fi
done

