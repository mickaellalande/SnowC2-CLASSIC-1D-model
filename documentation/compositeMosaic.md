# Tiling {#tiles}

1. @subpage compvsmosaic
2. @subpage dynTile


# Composite vs. Mosaic representations {#compvsmosaic}

CLASSIC can be run in two different configurations: composite or mosaic (see figure below)

In the composite mode, the structural vegetation attributes (including leaf area index, vegetation height, rooting depth) of PFTs that exist in a grid cell are averaged in proportion to their fractional coverages and then used in the grid-averaged energy and water balance calculations. As a result the entire grid cell is characterized by land surface physical environment (including soil temperature, soil moisture, fractional snow cover, and net radiation) that is common to all PFTs.

In contrast, in the mosaic mode, a grid box is split into multiple tiles representing individual PFTs (or land unit, such as soil texture e.g. Melton et al. (2017) \cite Melton2017-gp) for each of which energy, water and carbon balance calculations are performed separately. In principle, the mosaic mode may be used to represent tiles that are characterized by any chosen distinction such a lowlands vs. uplands, soil texture, vegetation, soil depth, etc.

As a result of the differences between the mosaic and composite configurations, the simulated carbon balance evolves somewhat differently in the two configurations despite being driven with identical climate forcing (see Melton and Arora (2014) \cite Melton2014-xk). For previous publications using this model capability see: Li and Arora (2012) \cite Li2012-f7f, Melton and Arora (2014) \cite Melton2014-xk, Shrestha et al. (2016) \cite Shrestha2016-do, and Melton et al. (2017) \cite Melton2017-gp.

Although pretty clever and powerful, running the mosaic version and interpreting the model results can be a logical nightmare, especially with competition on where the fractions of different tiles/mosaic change with time. So if you are new at this, please consider running the model in the composite mode.

Some initialization variables that impact upon composite vs. mosaic model runs:

- **FAREROT** Fractional coverage of mosaic tile on the modelled area
  - If you run with multiple tiles you need to specify what fraction of the grid cell the tile occupies.
-  MIDROT (<b>OBSOLETE - not used by CLASSIC v. 1.5+</b>) - Mosaic tile type identifier (1 for land surface, 0 for inland lake)
  - CLASSIC runs in the coupled models with a sub-grid lakes scheme (The Canadian Small Lake Model (CSLM); Verseghy and MacKay, 2017 \cite Verseghy2017-ys). The lakes can then be represented as a tile.


\image html "compVsMosaic_MeltonArora_BG_2014.png" "Schematic representation of the composite and mosaic approaches for the coupling of CLASS and CTEM models in a stand-alone mode. (From Melton and Arora, 2014)"
\image latex "compVsMosaic_MeltonArora_BG_2014.png" "Schematic representation of the composite and mosaic approaches for the coupling of CLASS and CTEM models in a stand-alone mode. (From Melton and Arora, 2014)"


# Dynamic tiling capability {#dynTile}

With CLASSIC version 1.5, the model has switched from dynamic to static allocation of its data structures (variables in src/ctemStateVars.f90 and src/classStateVars.f90). This change allows the model to have a variable number of tiles in use during a simulation. While the static allocation requires a set number of tiles at the start of the simulation (e.g. 20), the actual tiles that are being utilized can vary provided some code capabilities are added. Static allocation also results in faster model run times. A downside of this choice to move to static allocation is that you must specify `ignd`, `ican`, `icc`, and `l2max`, which were previously discovered from the model initial conditions file. Importantly - you must do a `make clean` to remove the model binary exectuable and objectfiles prior to compiling with new values for these parameters. A simple `make` won't update these since they are parameters in a module. The parameter values specified must match both your model initial conditions file and the chosen run parameters file.

These are the parameters to change (if needed) in src/classicParams.f90:

>  ! Originally read from netcdf initialization now made static\n
>  integer, parameter :: nlat = 1        !< Number of cells we are running, set in the code. Offline this always 1.\n
>  integer, parameter :: nmos = 1        !< Number of mosaic tiles, set in the code\n
>  integer, parameter :: ilg = nlat * nmos         !< nlat x nmos\n
>  integer, parameter :: ignd = 20        !< Number of soil layers, set in the code\n
>  \n
>  ! Originally read in from the run parameters file now made static\n
>  integer, parameter :: ican = 5             !< Number of CLASS (physics) pfts, set in the code\n
>  integer, parameter :: icc = 12              !< Number of CTEM (biogeochemical) pfts, set in the code\n
>  integer, parameter :: l2max = 3             !< Maximum number of level 2 CTEM PFTs, set in the code\n
