import yaml
import glob
import re
import os

# This will search for any string in the yaml based upon the key.
# E.g.
# key = 'biome' (must match those in the yamls; case insensitive)
# searchString = 'DBF' (case insensitive and will find all instances within the field)
# Returns a list of sites that match the searchString for that key.

# JM Jan 2023

cbcpath='/home/rjm001/CLASSIC/inputFiles/FLUXNETsites/'
#cbcpath='/home/rjm001/CLASSIC/inputFiles/FLUXNETsites_in_progress/COHERENT_Cinit/'
key = 'PFTs'
searchString = 'Sedge'  

for name in glob.glob(cbcpath + '*'):
    #print (name)
    site=name.rsplit('/', 1)[1]

    # Check if the site has an obs file.
    # if (not os.path.exists('/space/hall5/sitestore/eccc/crd/ccrp/users/scrd530/benchmark_classic/CBC/sites/daily_observations/FLUXNET2015_daily_obs_all_sites/'+site+'_obs.csv')):
    #     print ('Site:',site,'missing obs')
    #     if (os.path.exists('/space/hall5/sitestore/eccc/crd/ccrp/users/rjm001/fluxnet_sites/CLASSICv1.5tests/Ameriflux_beta_daily_obs_all_sites/'+site+'_obs.csv')):
    #         print ('Site:',site,'found in Ameriflux!')

    # Now check each site for the info we want from the yaml
    yamlfile = cbcpath + site + "/siteinfo.yaml"
    #try:
    with open(yamlfile, "r") as f:
        contents = yaml.safe_load(f)
    #print ("site in CBC", site)
    if (re.search(searchString,contents[key], re.IGNORECASE)):
        print (site," is ",contents[key])
