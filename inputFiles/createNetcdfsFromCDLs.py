import os, fnmatch
from os.path import splitext
import sys

# NOTE: This assumes you have nco installed and available (uses ncgen)

# JM Feb 2021

# This script is used to convert the cdl files tracked by git into 
# netcdf files for model ingestion. 

# This is the path where the cdl files are located. It will walk down
# the directory tree looking for the cdl files and then place the netcdf 
# in the same directory as the cdl file with the same naming convention.

# You could always do this manually, but this is convenient to do all sites at once.
dir_path = str(sys.argv[1])

def find_files(directory, pattern):
  for root, dirs, files in os.walk(directory):
    for basename in files:
      if fnmatch.fnmatch(basename, pattern):
        filename = os.path.join(root, basename)
        yield filename

for filename in find_files(dir_path, '??-???_init.cdl'):
  NCfilename = splitext(filename)[0] + '.nc'
  print ('creating ',NCfilename)
  cmd = "ncgen "+ filename + " -o " + NCfilename
  returned_value = os.system(cmd)
