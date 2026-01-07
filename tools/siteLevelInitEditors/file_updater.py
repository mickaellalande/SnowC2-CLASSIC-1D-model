#!/usr/bin/env python3

""" A general use text replacer. The work flow is as follows:
    -  in the same directory as this script, create a .txt file that contains the block of text you want to replace 
       and a .txt file containing the block of text that will be inserted.
    -  in the main method: (further comments there explaining what to do)
       -  specify the path or paths to the root directories containing the files (or subdirectories containing the files) that need to be updated.
       -  provide the names of your .txt files where appropriate.
       -  provide the file name that will be temporarily used for the updated files before the old versions are deleted.
       -  provide the string that will be used to find the names of the files that will be updated in the paths passed.
    -  the steps in the main method can be repeated as necessary."""

import os
import re
import xarray as xr

def directory_looper(paths, to_match, to_replace, to_add, updated_file, soilpH):
    """Loop through the paths passed, and update files with names that match contents of to_match."""
    for folder in paths:
        for root, dirs, files in os.walk(folder, topdown = False):
            for name in files:
                # Match file names.
                if re.match(to_match, name):
                    matched_file = os.path.join(root, name)
                    print(matched_file)
                    updated_file = os.path.join(root, updated_file)
                    try:
                        os.remove(updated_file)
                    except:
                        pass
                    print(updated_file)

                    # Find the soil pH
                    add_soilpH(matched_file, soilpH)
                    # Update the in.txt file
                    with open('in.txt', 'r') as f:
                        to_add = f.read()

                    # Create new file with the name contained in updated_file.
                    text_replacer(matched_file, to_replace, to_add, updated_file)

                    # Delete old matched file and rename the new file to the old file name.
                    file_replacer(matched_file, updated_file)
                    updated_file = os.path.basename(updated_file)
                    

def text_replacer(file, to_replace, to_add, updated_file):
    """Read original file, and create a new file that contains the updated contents of the original."""
    # Open a file to read.
    with open(file, 'r') as f:
        added = f.read()
        # Look for the text contained in to_replace and replace it with the text contained in to_add.
        added = added.replace(to_replace, to_add,1)
        # Create a new file with the name contained in updated_file and write all contents of the original file to it, but now with the updated text.
        with open(updated_file, 'w') as fw:
            fw.write(added)
    
def add_soilpH(file,soilpH):
    """Use the lat/lon to find the soil pH then add to file"""
    with open(file, 'r') as f:
        for line in f:
            if line.startswith("lat = "):
                latv = line.split(" ")[2]
            if line.startswith("lon = "):
                lonv = line.split(" ")[2]
    
        spH = soilpH['soilPH'].sel(lat=[latv], method='nearest').sel(
            lon=[lonv], method='nearest').sel(level=[0]).values[0][0][0]
        
        replac = "soilpH = " + str(spH) + ";"
        with open('in.txt', 'w') as fw:
            fw.write(replac)

def file_replacer(old_file, new_file):
    """Replace one file with another."""
    # Remove the file with the name contained in old_file.
    os.remove(old_file)
    # Rename the file with the name contained in new_file to have the name contained in old_file.
    os.rename(new_file, old_file)

        
def main():
    # Paths to files to be updated:
    paths = ["/home/rjm001/code/CLASSIC/inputFiles/FLUXNETsites"]
    # Initialize text vars.
    to_replace = ""
    to_add = ""

    #open soil pH file
    soilpH = xr.open_dataset(
        '/space/hall3/sitestore/eccc/crd/ccrp/scrd530/classic_inputs/initfiles/global/1deg/soil/OpenLandMap_soilPH_1deg.nc')
    
    # File containing text to be replaced:
    with open('out.txt', 'r') as f:
        to_replace = f.read()
    # File containing text to be inserted:
    with open('in.txt', 'r') as f:
        to_add = f.read()
    # Temporary name of the updated file(s):
    updated_file = "updated.cdl"
    # String on which to match the files that are being updated:
    to_match = ".*init\.cdl"
    # Loop through paths passed.
    directory_looper(paths, to_match, to_replace, to_add, updated_file, soilpH)
    
    
    # Repeat above steps as necessary:
    # with open('out2.txt', 'r') as f:
    #     to_replace = f.read()
    # with open('in2.txt', 'r') as f:
    #     to_add = f.read()
    
    #directory_looper(paths, to_match, to_replace, to_add, updated_file)
    


    
    
if __name__ == '__main__':
    main()
