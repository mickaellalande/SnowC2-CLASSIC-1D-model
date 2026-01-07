# Overview of CLASSIC {#overviewCLASSIC}

1. @subpage modelVers
2. @subpage modelstruct
2. @subpage PFTsCLASSIC
3. @subpage overviewCLASS "History of CLASS"
  - @subpage devHistory
4. @subpage overviewCTEM "History of CTEM"

# Model versions {#modelVers}

CLASSIC version 1.0 was released on October 31st, 2019 ([Zenodo repository](https://zenodo.org/record/3522407#.YrZOB3ZKiUk)) and described in two papers: Melton et al. 2020 \cite Melton2020-ss and Seiler et al. 2021 \cite Seiler2021-wq. From Melton et al.: 

> We  have  transitioned  CLASS-CTEM  from  an  internally developed model to an open-source community model, which we call the Canadian Land Surface Scheme including Biogeochemical  Cycles  (CLASSIC)  v.1.0.   CLASSIC  contains  many  technical  features  specifically  designed  to  encourage community use including software containerization for serial and parallel simulations, extensive benchmarking software and data (Automated Model Benchmarking; AMBER), self-documenting code, community standard formats for model inputs and outputs, amongst others. Here, we evaluate and benchmark CLASSIC against 31 FLUXNET sites where  the  model  has  been  tailored  to  the  site-level  condi- tions and driven with observed meteorology. Future versions of CLASSIC will be developed using AMBER and these initial benchmark results to evaluate model performance over time. CLASSIC remains under active development and the code, site-level benchmarking data, software container, and AMBER are freely available for community use.


CLASSIC version 1.5 (released DATE) introduced:
 * Non-structural carbohydrates (Asaadi et al., 2018 \cite Asaadi2018-xj)
 * An interactive nitrogen cycle (Asaadi and Arora, 2021 \cite Asaadi2021-yn)
 * An improved representation of biological nitrogen fixation (Kou-Giesbrecht and Arora, 2022 \cite Kou-Giesbrecht2022-wm)
 * Shrubs as a plant functional type (Meyer et al., 2021 \cite Meyer2021-xi)
 * Changed soil moisture interpolation scheme to Lagrange 1st order interpolation (Mackay et al. 2022) CITATION
 * Boreal peatland module fully integrated and operational (Wu et al. 2016 \cite Wu2016-zt)
 * Multilayer soil carbon scheme including turbation. Heterotrophic respiration routines rewrote to accommodate. 
 * Diffuse shortwave radiation as an additional model meteorological forcing
   * Can be used with two-leaf (shaded/sunlit) photosynthesis scheme
   * Can be used with four-band albedo scheme (In development - not recommended to use).
   * Added in read of four-band albedo lookup table
   * Added in read-in and use of black carbon deposition fluxes
 *  Benchmarking suite structure changed.
   * Now included in the repo are each site information file (`*.yaml`) and the CDL-format version of the model initial conditions file (`*_init.cdl`). There are also helper scripts to convert from the `*.cdl` format to `*.nc` and back (including nicer formatting than default from `ncdump`).
   * Number of sites greatly expanded.

Other code changes are listed in the [release notes](https://gitlab.com/cccma/classic/-/releases).


# Model structure {#modelstruct}

## main and mainCore

## Gather/Scatter

## Crops {#crops}

## Peatlands {#peatlandC}



# Plant Functional Types (PFTs) in CLASSIC {#PFTsCLASSIC}

The original PFT scheme of CLASSIC is as follows (parameters for this configuration are found in the configurationFiles/template_run_parameters.txt file):

| CLASS PFTs | CTEM PFTs --- | ---| ---|
|:---------------:|:---------:|----------------|------------|
| Needleleaf tree (NdlTr) | Evergreen ('NdlEvgTr') | Deciduous ('NdlDcdTr') |  |
| Broadleaf tree (BdlTr) | Evergreen ('BdlEvgTr') | Cold Deciduous ('BdlDCoTr') | Drought/Dry Decidous ('BdlDDrTr') |
| Crop (Crops) | \f$C_3\f$ ('CropC3  ') | \f$C_4\f$ ('CropC4  ')|  |
| Grass (Grass) | \f$C_3\f$ ('GrassC3 ') | \f$C_4\f$ ('GrassC4 ')|  |

Shrubs have been implemented as a fifth CLASS PFT allowing for the model physics to distinguish shrubs from trees. Parameters for this configuration are found in the configurationFiles/template_run_parameters_shrubs.txt file):

| CLASS PFTs | CTEM PFTs --- | ---| ---|
|:---------------:|:---------:|----------------|------------|
| Needleleaf tree (NdlTr) | Evergreen ('NdlEvgTr') | Deciduous ('NdlDcdTr') |  |
| Broadleaf tree (BdlTr) | Evergreen ('BdlEvgTr') | Cold Deciduous ('BdlDCoTr') | Drought/Dry Decidous ('BdlDDrTr') |
| Crop (Crops) | \f$C_3\f$ ('CropC3  ') | \f$C_4\f$ ('CropC4  ')|  |
| Grass (Grass) | \f$C_3\f$ ('GrassC3 ') | \f$C_4\f$ ('GrassC4 ')| Sedges ('Sedge   ')  |
| Broadleaf shrub (BdlSh) | Evergreen Shrubs ('BdlEvgSh') | Deciduous Shrubs ('BdlDCoSh')|  |

Canada specific PFTs have been implemented in CLASSIC to allow for a better representation of the carbon cycle within the canada domain (Curasi et al., 2022) \cite Curasi2022-ss. Parameters for this configuration are found in the configurationFiles/template_run_parameters_canada.txt file):

| CLASS PFTs | CTEM PFTs --- | ---| ---|
|:---------------:|:---------:|----------------|------------|
| Needleleaf tree (NdlTr) | Evergreen ('NdlEvgTr') | Deciduous ('NdlDcdTr') | Continental Evergreen ('NdlEvgCtTr') | Interior Evergreen ('NdlEvgItTr') |
| Broadleaf tree (BdlTr) | Evergreen ('BdlEvgTr') | Cold Deciduous ('BdlDCoTr') | Drought/Dry Decidous ('BdlDDrTr') |
| Crop (Crops) | \f$C_3\f$ ('CropC3  ') | \f$C_4\f$ ('CropC4  ')|  |
| Grass (Grass) | \f$C_3\f$ ('GrassC3 ') | \f$C_4\f$ ('GrassC4 ')| Sedges ('Sedge   ')  |
| Broadleaf shrub (BdlSh) | Evergreen Shrubs ('BdlEvgSh') | Deciduous Shrubs ('BdlDCoSh')|  |

The model is not limited to these PFTs provided PFT-specific parameters can be determined. As well it is important to thoughtfully integrate new PFTs. Within the model code there are checks for unknown PFTs based upon the short names of the known PFTs in the tables above.

# Rate change equations for carbon pools {#CTEMRateChgEqns}


From the gross canopy photosynthesis rate (\f$G_{\text{canopy}}\f$, @ref PHTSYN3.f), maintenance and growth respirations (\f$R_\mathrm{m}\f$ and \f$R_\mathrm{g}\f$, @ref mainres.f), and
heterotrophic respiration components (\f$R_{\text{h,H}}\f$ and \f$R_{\text{h,D}}\f$, @ref hetres_mod.f90), it is possible to estimate the change in carbon amount of the model's five pools.

When the daily NPP (\f$ G_{canopy} - R_\mathrm{m} - R_\mathrm{g}\f$) is positive, carbon is allocated to the plant's live carbon pools and the rate of change is given by

  \f[
  \frac{\mathrm{d}C_i}{\mathrm{d}t} = a_{fi} \left(G_{canopy}-R_\mathrm{m}-R_\mathrm{g} \right) - D_i - H_i - M_i \\ \quad i = {L, S, R} \qquad (Eqn 1)
   \f]
   <!-- {#rate_change_eqns_live_pools} -->

where \f$a_{fi}\f$ is the corresponding allocation fractions for each pool (stem, root and leaves) and \f$D_i\f$ is the litter produced from these components as explained in @ref phenolgy.f90. \f$H_i\f$ is the loss associated with fire that releases \f$CO_2\f$ and other trace gases to the atmosphere and \f$M_i\f$ is the mortality associated with fire that contributes to the litter pool as explained in @ref disturb.f90.

If the daily NPP is negative (\f$G_{canopy} < R_\mathrm{m}\f$, \f$R_\mathrm{g} = 0\f$), the rate of change is given by

\f[
 \frac{\mathrm{d}C_i}{\mathrm{d}t} = a_{fi}G_{canopy} - R_{m,i}  - D_i  - H_i - M_i, \\ \quad i = {L, S, R} \qquad (Eqn 2)
 \f]
<!-- \label{rate_change_eqns_live_pools2} -->

Negative NPP causes the plant to lose carbon from its live carbon pools due to respiratory costs in addition to the losses due to litter production (\f$D_i\f$) and disturbance (\f$H_i\f$, \f$M_i\f$).

The rate change equations for the litter and soil carbon pools are given by
\f{eqnarray*}{
\frac{\mathrm{d}C_\mathrm{D}}{\mathrm{d}t} &=& D_\mathrm{L} + D_\mathrm{S} +
D_\mathrm{R} + M_\mathrm{L} + M_\mathrm{R} + M_\mathrm{S} - H_\mathrm{D} -C_{\mathrm{D} \rightarrow \mathrm{H}} - R_{h,D} \\
\frac{\mathrm{d}C_\mathrm{H}}{\mathrm{d}t} &=& C_{\mathrm{D} \rightarrow
\mathrm{H}} - R_{h,H}
\qquad (Eqn 3)
\f}
<!-- \label{rate_change_eqns_dead_pools}, -->

where \f$C_{\mathrm{D} \rightarrow \mathrm{H}}\f$ represents the transfer of
humified litter to the soil carbon pool and \f$H_\mathrm{D}\f$
is loss associated with burning of litter associated with fire that releases
\f$CO_2\f$ and other trace gases to the atmosphere.


# Historical overview of the Canadian Land Surface Scheme (CLASS) {#overviewCLASS}

The Canadian Land Surface Scheme, CLASS, was originally developed for use with the Canadian Global Climate Model (CanGCM) (Verseghy, 1991 \cite Verseghy1991-635 ; Verseghy et al., 1993 \cite Verseghy1993-1ee ). The table at the end of this overview summarizes the development of CLASS from the late 1980’s onward.

The basic function of CLASS is to integrate the energy and water balances of the land surface forward in time from an initial starting point, making use of atmospheric forcing data to drive the simulation. When CLASS is run in coupled mode with a global or regional atmospheric model, the required forcing data are passed to it at each time step over each modeled grid cell from the atmospheric driver. CLASS then performs its internal calculations, evaluating a suite of prognostic and diagnostic variables such as albedo and surface radiative and turbulent fluxes, which are in turn passed back to the driver. CLASS can also be run in uncoupled or offline mode, using forcing data derived from field measurements, and the output values of its prognostic and diagnostic variables can then be validated against observations.

CLASS models separately the energy and water balances of the soil, snow, and vegetation canopy (see the diagram below). The basic prognostic variables consist of the temperatures and the liquid and frozen moisture contents of the soil layers; the mass, temperature, density, albedo and liquid water content of the snow pack; the temperature of the vegetation canopy and the mass of intercepted rain and snow present on it; the temperature and depth of ponded water on the soil surface; and an empirical vegetation growth index (which is not used when the biogeochemistry module, CTEM, is turned on). These variables must be initialized, and a set of physical parameters describing the soil and vegetation existing on the modelled area must be assigned background values, at the beginning of the simulation (see @ref forcingData).

At each time step, CLASS calculates the bulk characteristics of the vegetation canopy on the basis of the vegetation types present over the modelled area. In a pre-processing step, each vegetation type is assigned representative values of parameters such as albedo, roughness length, annual maximum and minimum plant area index, rooting depth and so on (see @ref initProgVar). These values are then aggregated over vegetation categories identified by CLASS (commonly includes: needleleaf trees, broadleaf trees, shrubs, crops, and grass, i.e. short vegetation). The physiological characteristics of the vegetation in each category are determined at the current time step using the aggregated background parameters and assumed annual or diurnal variation functions (unless provided by CTEM). These physiological characteristics are then aggregated to produce the bulk canopy characteristics for the current time step.

\image html "schematicDiagramOfClass.png" "Schematic Diagram Of CLASS"
\image latex "schematicDiagramOfClass.png" "Schematic Diagram Of CLASS"

In performing the surface flux calculations the modeled area is divided into up to four subareas: bare soil, vegetation over soil, snow over bare soil, and vegetation over snow. The fluxes of these sub-regions are determined each CLASS timestep, the average value over the modelled area is found and the sub-region fluxes are initialized with the average value at the start of the next time step. The fractional snow coverage is determined using the concept of a threshold snow depth. If the calculated snow depth is less than this value, the snow depth is set to the threshold value and the fractional snow cover is calculated on the basis of conservation of snow mass, i.e. the fractional cover decreases since the snow depth is set to the threshold value. The fluxes are calculated for each of the four subareas, and these and the prognostic variables are then areally averaged before being passed back to the atmospheric model.

Originally CLASS performed only one set of these calculations for each grid cell of the model domain. In more recent versions, a “mosaic” option has been added to handle sub-grid scale heterogeneity more effectively by allowing tiling. When this option is utilized, each grid cell is divided into a user-specified number of “tiles”, and the CLASS calculations are performed in turn over each. The surface fluxes are averaged, but the prognostic variables are kept separate for each of the tiles of the mosaic between time steps (see [here for more](@ref compvsmosaic)).

In the CLASSIC offline driver, a gather-scatter operation is included in the driver, mimicking the practice in atmospheric models of “gathering” land surface points on latitude circles onto long vectors prior to the calculations (e.g. src/ctemGatherScatter.f90 or src/classGatherScatter.f90), for improved computational efficiency on vector supercomputers. For CLASS, the mosaic tiles on each of the modelled grid cells are “gathered” onto long arrays prior to calling the CLASS subroutines (thus collapsing the first two dimensions of the arrays into one), and subsequently “scattered” back onto the grid cells before performing the diagnostic averaging calculations. Future code developments will work to make this necessity easier to work with, however developments will follow those adopted in the CanESM framework.

## Development timeline of CLASS {#devHistory}

\f[
\begin{array}{ | c | c | l | }
1.0 & \text{April 1989} & \text{Basic thermal and hydrological model of snow and soil.} \\
2.0 & \text{August 1991} & \text{Addition of vegetation thermal and hydrological model.} \\
2.1 & \text{May 1993} & \text{Full vectorization of code to enable efficienr running on vector supercomputers.} \\
2.2 & \text{April 1994} & \text{Augmentation of diagnostic calculations; incorporation of in-line comments throughout;} \\
    & & \text{development of a parallel stand-alone version of the model for use with field data.} \\
2.3 & \text{December 1994} & \text{Revisions to diagnostic calculations; new near-surface atmospheric stability functions.} \\
2.4 & \text{August 1995} & \text{Complete set of water budget diagnostic calculations; parametrizations of organic soils} \\
    & & \text{and rock soils; allowance for inhomegeneity between soil layers; incorporation of variable} \\
    & & \text{surface detention capacity.} \\
2.5 & \text{January 1996} & \text{Completion of energy budget diagnostic calculations.} \\
2.6 & \text{August 1997} & \text{Revisions to surface stability function calculations.} \\
2.7 & \text{December 1997} & \text{Incorporation of variable soil permeable depth; calculation of soil thermal and hydraulic} \\
    & & \text{properties based on textural composition; modified surface temperature iteration scheme.} \\
3.0 & \text{December 2002} & \text{Improved treatment of soil evaporation; complete treatment of organic soils; new canopy} \\
    & & \text{conductance formulation; preliminary routines for lateral movement of soil water; enhanced} \\
    & & \text{snow density and snow interception; improved turbulent transfer from vegetation; mosaic} \\
    & & \text{formulation.} \\
3.1 & \text{April 2005} & \text{Faster surface temperature iteration scheme; refinements to leaf boundary resistance} \\
    & & \text{formulation; improved treatment of snow sublimation and interception; transition to} \\
    & & \text{Fortran 90 and single precision variables.} \\
3.2 & \text{May 2006} & \text{Option for multiple soil layers at depth; additional liquid water content of snow pack;} \\
    & & \text{revised radiation transmission in vegetation.} \\
3.3 & \text{December 2006} & \text{Separate temperature profile curve fit for snow and soil; multiple-layer option for ice} \\
    & & \text{sheets; water and energy balance checks for each time step; modifications to soil hydraulic} \\
    & & \text{conductivity calculations.} \\
3.4 & \text{April 2008} & \text{Streamline and clean up code; updated soil thermal conductivity calculations; revisions to} \\
    & & \text{handling of water stored on vegetation.} \\
3.5 & \text{December 2010} & \text{Updated field capacity calculation; revised treatment of water on canopy; reworked} \\
    & & \text{calculation of baseflow.} \\
3.6 & \text{December 2011} & \text{Revised ponding depth over organic soils; revised snow albedo refreshment threshold; new} \\
    & & \text{snow thermal conductivity algorithm; interface with Canadian Terrestrial Ecosystem Model} \\
    & & \text{(CTEM).} \\
3.6.1 & \text{December 2016} & \text{New treatment of bare soil albedo; new optional four-band snow albedo formulation;} \\
    & & \text{fixes to guard against overshoots in water drawdown by evapotranspiration; upper limit on} \\
    & & \text{snow depth.} \\
3.6.2 & \text{July 2019} & \text{CLASS is formerly incorporated into the Canadian Land Surface Scheme including} \\
    & & \text{Biogeochemical Cycles (CLASSIC).} \\
\end{array}
\f]


# Historical overview of the Canadian Terrestrial Ecosystem Model (CTEM) {#overviewCTEM}

Version 1 of the CTEM is the terrestrial carbon cycle component of the second generation Canadian Earth System Model (CanESM2) (Arora et al., 2011)\cite Arora2011-79f where it was coupled to version 2.7 of the Canadian Land Surface Scheme (CLASS). CTEM v. 2.0 (Melton and Arora, 2016) \cite Melton2016-zx has been coupled to CLASS v. 3.6 (Verseghy, 2012) \cite Verseghy2012-c0e. Together CLASS and CTEM form the CLASSIC model which is capable of being run online in the CanESM family of models or offline, driven by observation-based meteorological forcings. CLASSIC models terrestrial ecosystem processes by tracking the flow of carbon, and optionally nitrogen, through three living vegetation components (leaves, stem and roots) and two dead carbon pools (litter and soil).

<!-- \f[
\begin{table}[]
\caption{CTEM and peatland PFTs and their mapping to the CLASS PFTs}
\label{my-label}
\begin{array}{|l|l|l|l|l|l|}
\hline
\multicolumn{1}{|c|}{CLASS PFTs} & \multicolumn{3}{c|}{CTEM PFTs}                    & \multicolumn{2}{c|}{Peatland PFTs}  \\ \hline
Needleleaf tree                  & Evergreen & Deciduous      &                      &                  &                  \\ \hline
Broadleaf tree                   & Evergreen & Cold Deciduous & Drought/Dry Decidous & Evergreen Shrubs & Deciduous Shrubs \\ \hline
Crop                             & C$_3$     & C$_4$          &                      &                  &                  \\ \hline
Grass                            & C$_3$     & C$_4$          &                      & Sedges           &                  \\ \hline
\end{array}
\end{table}
\f] -->

The amount of carbon and nitrogen in these five carbon pools is simulated prognostically (see below). In the CLASSIC framework, CLASS uses structural vegetation attributes (including LAI, vegetation height, canopy mass and rooting depth) simulated by CTEM, and CTEM uses soil moisture, soil temperature and net radiation calculated by CLASS. Combined, CLASS and CTEM simulate the atmosphere--land fluxes of energy, water,  \f$CO_2\f$, nitrogen, and \f$CH_4\f$.

Version 1.0 of CTEM is described in a collection of papers detailing parametrization of photosynthesis, autotrophic and heterotrophic respiration (Arora, 2003) \cite Arora2003-3b7; phenology, carbon allocation, biomass turnover and conversion of biomass to structural attributes (Arora and Boer, 2005) \cite Arora2005-6b1; dynamic root distribution (Arora and Boer, 2003) \cite Arora2003838; and disturbance (fire) (Arora and Boer, 2005) \cite Arora20052ac. These processes are modelled over prescribed fractional coverage of (at the time, typically) nine PFTs (Wang et al., 2006) \cite Wang2006-he and determine the structural vegetation dynamics including vegetation biomass, LAI, vegetation height, fraction of roots in each of the three soil layers, leaf onset and offset times and primary \f$CO_2\f$ fluxes of gross primary productivity (GPP) and NPP.

CTEM v. 2.0 is described in Melton and Arora, 2016 \cite Melton2016-zx. CTEM v. 2.0 can be run in two different modes, either (i) using specified fractional coverage of its PFTs, or (ii) allowing the fractional coverage of its non-crop PFTs to be dynamically determined based on competition between PFTs. The parametrization for simulating competition between PFTs is found in @ref competition_mod.f90 and described in Arora and Boer (2006a,b), \cite Arora2006-ax \cite Arora2006-pp,  Melton and Arora, 2016 \cite Melton2016-zx, and Shrestha et al. (2016) \cite Shrestha2016-do. The fire parametrization has also been refined in the new model version as described in @ref disturb.f90 and Melton and Arora, 2016 \cite Melton2016-zx.