# How do I ... / Something has gone wrong! {#howDoI}

1. @ref chgSoil
2. @ref failStart
3. @ref makeClean
4. @ref soilColourIndex
5. @ref runInWindows
6. @ref windowsEndings
7. @ref macrun
8. @ref izrefsite
9. @ref reff0

----

# Change the number/depth/etc of soil layers? {#chgSoil}

Information about the soil layers is taken in from the initialization netcdf file. model_state_drivers::read_modelsetup reads the netcdf for the number of soil layers (*ignd*). The *DELZ* variable in the intialization file is the thickness of each layer. So change the initialization file for the number and thicknesses of the soil layers as required and CLASSIC will simply read them in.

# My run starts and I get an immediate fail! {#failStart}

If you start a new run and it immediately fails with an output something like:

        Singularity jormelton-containerCLASSIC-master-latest.simg:~/Documents/CLASSIC> bin/CLASSIC configurationFiles/template_job_options_file.txt 254.85/53.98
         in process           0           1           1   254.84999999999999        53.979999999999997                1           1           1           1
          CANOPY ENERGY BALANCE         1       8092.77620036       8077.49643236
                 0.000000       1.875648      -2.045030      -0.189964       0.000000    8073.385790
                 0.000000       1.339229     257.925000
         died on   254.84999999999999        53.979999999999997 
         
This is a general indication of a problem in your initialization file. Often it is fixed by setting snow/liquid in the canopy to zero (see @ref initPhysProgVar). You may also get snow energy/water balance failing, again look to @ref initPhysProgVar for advice with setup.

# My run won't compile but I didn't do anything I can think of. {#makeClean}

Sometimes you need to run a '`make mode=???? clean`' (where mode is the one you used to compile the code, for site-level it is typically `serial`) to clean out old .mod and .o files that might be causing issues. Do that and then recompile.

# initFileConverter tells me I need a soil colour index {#soilColourIndex}

Soil colour index is described in @ref soilData and used in soilProperties.f90. A global file of soil colour index produced by Peter Lawrence (NCAR) can be obtained from ftp://ftp.cccma.ec.gc.ca/pub/jmelton/mksrf_soilcol_global_c090324.nc

# Run on a Windows Machine {#runInWindows}

Two options here. The first uses a VirtualBox and the second uses the Windows Subsystem for Linux (WSL). Please choose your flavour depending on your circumstances. Either solution will not be quite as fast a running on a native Linux machine.

## VirtualBox

Courtesy of E. Humphreys. 

1. Install Oracle VirtualBox
2. Install Ubuntu as a virtual machine with ~ 100GB drive space
3. Install Dropbox on Ubuntu (to share files between machines but there are other ways to do this)
4. Install Singularity.  
  1. In Terminal, type: `sudo apt install singularity-container`
5. Download our container into a working directory that is shared via Dropbox with the windows machine, e.g. `~/Dropbox/CLASSIC_working/`
6. Shell into the container:  `singularity shell CLASSIC_container.simg`
7. In this container, get CLASSIC code and folders:  

**Note the text files need to have Unix line endings if created in Windows (see @ref windowsEndings), make sure to convert to Unix (LF) from Windows (CR LF).**

## Windows Subsystem for Linux (WSL) 

Courtesy of A. Mavrovic.
 
1. Install WSL 2 following Microsoft’s instructions (https://docs.microsoft.com/en-us/windows/wsl/install-win10) 
2. Be sure to update WSL 1 to WSL 2 as CLASSIC requires functionalities only available in WSL 2 
3. Install your Linux distribution of choice from Microsoft Store (Ubuntu available) 
4. (Optional) Install Windows Terminal from Microsoft Store. By default, Windows Terminal open the PowerShell but it can be set to Ubuntu (open the parameters file). You should now be able to open Linux anywhere in your Windows system. 
5. Install CLASSIC from your Linux terminal as if you were on a Linux machine. 
6. If you struggle to install Singularity, try the installation procedure proposed by Sylabs (https://sylabs.io/guides/3.0/user-guide/installation.html) 

Additional notes by S.R. Curasi.

1. Singularity doesn't work on WSL1 so you must use WSL2.
2. Running the "cow" example (https://singularity-tutorial.github.io/02-basic-usage/) will help to diagnose if singularity is working correctly.
3. Singularity and classics won't work 100% of the time if the model run files are located in the Windows file system as opposed to the Linux file system created by WSL2, which is basicually just a folder within windows. If you're at /mnt/c/Users/ in WSL2 that means you're in the Windows file system. Bringing files from the windows system over to the Linux system might also cause problems so it's better to pull them straight in through WSL2.
4. After I installed WSL2 I had a lot of persistent network connection issues within WSL2, which caused further issues with singularity and installing software. Some combination of ping www.google.com, sudo apt-get install --reinstall resolvconf and opening /etc/resolv.conf and editing nameserver to 8.8.8.8 resolved these issues.

# Gotcha for Windows users {#windowsEndings}

Text files created on DOS/Windows machines have different line endings than files created on Unix/Linux. DOS uses carriage return and line feed ("\r\n") as a line ending, which Unix uses just line feed ("\n"). You need to be careful about transferring files between Windows machines and Unix machines to make sure the line endings are translated properly.(from http://www.cs.toronto.edu/~krueger/csc209h/tut/line-endings.html)

Any CLASSIC tool that reads in ASCII expects the Linux/Unix line endings.

# Run on a mac via virtual machine {#macrun}

This was adapted from https://sylabs.io/guides/3.0/user-guide/installation.html by S.R. Curasi and tested on macOS Catalina 10.15.7

1. Install Homebrew, virtual box, and other dependencies
  1. in terminal, type: >/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
  2. in terminal, type: >brew install virtualbox --cask
  3. in terminal, type: >brew install vagrant --cask
  4. in terminal, type: >brew install vagrant-manager --cask
2. Install singularity in the virtual machine and start it
  1. in terminal, type: >mkdir vm-singularity
  2. in terminal, type: >cd vm-singularity
  3. in terminal, type: >vagrant destroy
  4. in terminal, type: >rm Vagrantfile
  5. in terminal, type: >export VM=sylabs/singularity-3.0-ubuntu-bionic64
  6. in terminal, type: >vagrant init $VM
  7. in terminal, type: >vagrant up
  8. in terminal, type: >vagrant ssh
3. check if singularity is installed
  1. in terminal, type: >singularity version

Note 1: To copy files to the VM you can use the example code below (the default password which is "vagrant"). It’s also possible to copy files using the shared vm-singularity folder, or pulling them directly into the VM from git.
  1. in terminal, type: >scp -P 2222 -o StrictHostKeyChecking=no -r CLASSIC_container.simg vagrant@127.0.0.1:

Note 2: It’s easy to set up vs code to work with the VM on mac using the instructions at https://medium.com/@lopezgand/connect-visual-studio-code-with-vagrant-in-your-local-machine-24903fb4a9de

# My site-level run fails with IZREF = 1 but is okay with IZREF = 2 {#izrefsite}

Site-level runs should use IZREF = 1 with appropriate values of ZRFM and ZRFH for your site (described in @ref setupJobOpts). However it is possible, especially with CTEM (biogeochemistry) on that the vegetation will grow taller than ZRFM/ZRFH causing the model to crash. The model may simulate vegetation at the site taller than reality due to prior land use history at the site, model bias, etc. In this case you will need to raise the ZRFM/ZRFH heights to remain above the canopy. 

# Strange compiler fail due to reff0 #reff0

If you are compiling and you get a fail similar to below:

        gfortran -O0 -g -fdefault-real-8 -ffree-line-length-none -finit-real=snan -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -JobjectFiles  -I/fs/ssm/hpco/exp/hdf5-netcdf4/serial/static/gnu-4.8.5/01/ubuntu-18.04-amd64-64/include -c src/fourBandAlbedo.f90 -o objectFiles/fourBandAlbedo.o
          src/fourBandAlbedo.f90:79:67:

                 GCGAT, REFF0_LAND, DELT, ILG, 1, NML, GCMIN, GCMAX)
                                                                   1
        Error: Type mismatch in argument ‘reff0’ at (1); passed REAL(4) to REAL(8)
        src/fourBandAlbedo.f90:82:78:

            call RDSWET(SNOGAT, WSNOGAT, REFGAT, GCGAT, REFF0_LAND, DELT, ILG, 1, NML)
                                                                              1
        Error: Type mismatch in argument ‘reff0’ at (1); passed REAL(4) to REAL(8)
        src/fourBandAlbedo.f90:86:77:

                 REFF0_LAND, DELT, ZSNMIN, ZSNMAX2, ILG, 1, NML, GCMIN, GCMAX)
                                                                             1
        Error: Type mismatch in argument ‘reff0’ at (1); passed REAL(4) to REAL(8)
        Makefile:199: recipe for target 'objectFiles/fourBandAlbedo.o' failed
        make: *** [objectFiles/fourBandAlbedo.o] Error 1

It could be due to VSCode creating a `.mod` file in the `/src` directory. Remove any `.mod` files in `src` and recompile then the problem should resolve. It appears to be a bug in VSCode.