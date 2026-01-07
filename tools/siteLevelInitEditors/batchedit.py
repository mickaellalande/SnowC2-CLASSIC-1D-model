import os, fnmatch
from os.path import splitext
import sys

dir_path = str(sys.argv[1])

def find_files(directory, pattern):
  for root, dirs, files in os.walk(directory):
    for basename in files:
      if fnmatch.fnmatch(basename, pattern):
        filename = os.path.join(root, basename)
        yield filename

# for filename in find_files(dir_path, 'siteinfo.yaml'):
#   cmd = "sed -i '/^notes:.*/i PFTs: \"\" # use same PFT naming convention as in run_parameters.txt/model code.\\nleap: \".false.\" #.true. if the meteorology contains leap years'"+" "+ filename
#   print (cmd)
#   returned_value = os.system(cmd)

for filename in find_files(dir_path, '??-???_init.cdl'):
  cmd = "sed -i '/ALBS:long_name/c\\\t\tALBS:long_name = \"Snow albedo\" ;'"+" "+ filename
  print (cmd)
  returned_value = os.system(cmd)

