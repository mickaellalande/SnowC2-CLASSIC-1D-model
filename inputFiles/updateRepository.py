import os, fnmatch
from os.path import splitext
import shutil
import subprocess
import sys

# NOTE: This calls a script which assumes you have nco installed 
# and available (uses ncdump) along with the python library xarray!

# This script is used to ensure we are able to keep track of changes in the 
# benchmark files. We don't want to track the netcdf files themselves, but 
# it is very easy to convert them to and from cdl format. We also use a 
# script that makes the cdl's nicely formatted for manual editing. This 
# script here prompts the user to make sure they have done things 
# in the right order then:
# 1. generates the cdl files of the site-level inits
# 2. does a dry run of the commit so the user can be sure of what is going 
# to go in the commit
# 3. then performs the commit after prompting for the commit message
# Set the upper directory for the flux sites.
dir_path = str(sys.argv[1])

# Point to the nc_to_cdl.py script that is included in
# the CLASSIC repo, if you haven't moved it, this should work.
nc_to_cdl="../tools/siteLevelInitEditors/nc_to_cdl.py"

passon = False

print('The following questions are to make sure no mistakes were made.')
changelog=input("Did you edit the CHANGELOG.md for your changes? (y/n)")
if changelog == "y":
  newfilename = input("Does the newly edited file retain the old name structure (??-???_init.nc) (n/y)")
  if newfilename == "y":
    safetocdl = input("Is it safe to generate new cdl files from the init netcdfs, i.e. you don't have any edited cdls lying around?(n/y)")
    if safetocdl == "y":
      passon=True

def find_files(directory, pattern):
  for root, dirs, files in os.walk(directory):
    for basename in files:
      if fnmatch.fnmatch(basename, pattern):
        filename = root
        yield filename

if passon:
  for filename in find_files(dir_path, '??-???_init.nc'):
    cmd = "python "+ nc_to_cdl +" "+ filename
    returned_value = os.system(cmd)
  print ('DRY RUN:=====')
  gitcmd = "git add . --dry-run"
  rungitv = os.system(gitcmd)
  
  seemok = input("Do you wish to commit these changes to the repo?(y/n)")
  if seemok == "y":
    gitcmd2 = "git add ."
    rungit2 = os.system(gitcmd2)
    gitmsg = input("Enter your git commit message:")
    rungit3 = subprocess.run(["git", "commit", "-a","-m",gitmsg])
  else:
    print('commit aborted')
    quit()
else:
    quit()
