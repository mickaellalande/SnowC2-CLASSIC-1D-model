# builtin libs
from os.path import basename, splitext, abspath
import os, re
import sys
from subprocess import check_output

# external libs
import xarray as xr

# NOTE: This assumes you have nco installed and available (uses ncdump)! 

indir = str(sys.argv[1])

# String on which to match the files that are being updated:
to_match = ".*init.nc"

DO_ROUNDING = False # swap this to True if you want rounding to occur
DECIMAL_PRECISION = 3

for root, dirs, files in os.walk(indir, topdown = True):
    for name in files:
        #print (name)
        # Match file names.
        if re.match(to_match, name):
            IN_NC = os.path.join(root, name)
            OUT_CDL = splitext(IN_NC)[0] + '.cdl'
            with xr.open_dataset(IN_NC) as ds:
                # pull the header from the original file
                nc_header = check_output('ncdump -h {nc}'.format(nc=IN_NC), shell=True).decode('utf-8')
                # reformat the variables/dimensions for easier editing
                with open(OUT_CDL, 'w') as f:
                    # it might not be necessary to correct the name at the top of the CDL but we do it anyway
                    f.writelines(nc_header.replace('}\n', '').replace(splitext(basename(IN_NC))[0], splitext(basename(OUT_CDL))[0]))
                    f.writelines('data:\n')
                    # not just data-vars but dimensions as well
                    for var in ds.variables:
                        # squeeze() gets rid of extraneous array dims
                        squeezed = ds[var].squeeze()
                        arr_data = squeezed.data
                        desc = squeezed.attrs.get('long_name') if 'long_name' in squeezed.attrs.keys() else squeezed.attrs.get('standard_name')
                        # 2D arrs: convert to comma-separated table-ish
                        if len(squeezed.shape) == 2:
                            f.writelines(f'// {desc} ({", ".join(squeezed.dims)})\n')
                            f.writelines(f'{var} = \n')
                            rows, cols = squeezed.shape
                            last_idx = (rows - 1) * (cols - 1)
                            for i in range(rows):
                                for j in range(cols):
                                    val = round(arr_data[i, j], DECIMAL_PRECISION) if DO_ROUNDING else arr_data[i, j]
                                    if i * j == last_idx:
                                        # last value must end in a semicolon
                                        f.write(f'{val} ;')
                                    else:
                                        f.write(f'{val}, ')
                                f.writelines('\n')
                        # 1D arrs: output as comma-separated one line with semicolon at the end
                        elif len(squeezed.shape) == 1:
                            f.writelines(f'// {desc} ({", ".join(squeezed.dims)})\n')
                            f.writelines(f'{var} = \n')
                            nums = squeezed.shape[0]
                            for i in range(nums):
                                val = round(arr_data[i], DECIMAL_PRECISION) if DO_ROUNDING else arr_data[i]
                                if i == nums - 1:
                                    f.write(f'{val} ;')
                                else:
                                    f.write(f'{val}, ')
                            f.writelines('\n')
                        # 0D arrs: output single-value with semi-colon at the end
                        elif len(squeezed.shape) == 0:
                            val = arr_data.round(DECIMAL_PRECISION) if DO_ROUNDING else arr_data
                            f.writelines(f'// {desc}\n')
                            f.writelines(f'{var} = {val} ;\n')
                        f.writelines('\n')
                    f.writelines('}\n')
            print(f'Created CDL file: {abspath(OUT_CDL)}')
