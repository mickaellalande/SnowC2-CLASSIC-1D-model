# The Canadian Land Surface Scheme Including biogeochemical Cycles (CLASSIC)

⚠️ **Important Notice**  
This repository is **not the official CLASSIC repository**.  
It is a **research copy** prepared specifically for the publication of the paper:

> *Improving the CLASSIC Snow Model to Better Simulate Arctic Snowpacks*  
> Submitted to *Geoscientific Model Development (GMD)*

The official CLASSIC v1.0 release is available here:  
👉 [https://gitlab.com/cccma/classic](https://gitlab.com/cccma/classic) (last access: 2025-09-02)  

For information about the development version of CLASSIC, please contact:  
**Joe Melton** – Joe.Melton@ec.gc.ca  

Further details on CLASSIC itself are available on the [CLASSIC webpage](https://cccma.gitlab.io/classic_pages/).


## Purpose of This Repository

This repository contains a modified version of CLASSIC that implements **six new model developments** aimed at improving snow simulations in Arctic conditions. These modifications were developed and tested as part of the study submitted to GMD.

The large files (including binaries, input, and output files) are furnished in a separate **Zenodo repository**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18157337.svg)](https://doi.org/10.5281/zenodo.18157337)


The analysis scripts (for data processing and figure generation) are provided separately in the repository: **[SnowC2-CLASSIC-1D-analysis](https://github.com/mickaellalande/SnowC2-CLASSIC-1D-analysis/tree/v0.1.0)**


## Branches in This Repository

- **[SnowC2-1D](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/tree/SnowC2-1D)**  
  Full implementation of all new developments, including job option files to launch the model.  

- **[SnowC2-1D-clean](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/tree/SnowC2-1D-clean)**  
  Cleaned version containing only the modified files.  
  Intended for direct [comparison](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/compare/develop-1.8...SnowC2-1D-clean) against the `develop-1.8` branch.  

- **[SnowC2-1D-backup](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/tree/SnowC2-1D-backup)**  
  Legacy version including the earlier compaction scheme (not factorized) and other tests.  

- **[develop-1.8](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/tree/develop-1.8)**  
  Starting point for this work, originating from a private GitLab repository maintained by Joe Melton.  


## Model Developments

### 1. Physics Developments

An overall (non-official) draft description of the snow model (before our model developments) is accessible in this Google Doc: https://docs.google.com/document/d/14TArlBpEcOGD1KmYKjZh603KAYY6V74ilX7Uv3zo2oY/edit?usp=sharing (last access: 2025-09-03)

The official documentation of CLASSIC v1.0 is accessible here: https://cccma.gitlab.io/classic/index.html (last access: 2025-09-03)

---
#### 1.1. Thermal conductivity at the top of the first soil layer  

Pass `ZERO` to [soilHeatFluxPrep.f90](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/soilHeatFluxPrep.f90) instead `ZSNOW` to the subfractions without snow in [energyBudgetDriver.f90](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/energyBudgetDriver.f90):
- *CANOPY OVER BARE GROUND*: https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blame/SnowC2-1D/src/energyBudgetDriver.f90#L1357
- *BARE GROUND*: https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blame/SnowC2-1D/src/energyBudgetDriver.f90#L1556

See more information at: https://gitlab.com/cccma/classic/-/issues/119 (last access: 2025-09-03)

---

#### 1.2. Bottom snow temperature 


In [snowTempUpdate.f90](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/snowTempUpdate.f90), the computation of the bottom snow temperature $T_{s, b}$ has been modified from:

$$
T_{s,b} = \frac{d_s T_s + d_1 T_1}{d_s + d_1}
$$

to

$$
T_{s,b} = 
\frac{\dfrac{T_s}{d_s} + \dfrac{T_1}{d_1}}
     {\dfrac{1}{d_s} + \dfrac{1}{d_1}},
$$

where $T_1$ is the first soil layer temperature, $T_s$ the snow temperature, and $d_1$ and $d_s$ their respective thickness.


- https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/snowTempUpdate.f90#L152

---

#### 1.3. Windless exchange coefficient  

Incorporate a windless transfer coefficient $E_0$ into the sensible heat flux $Q_H$ (W m $^{-2}$) calculation of 2 W m $^{-2}$ K $^{-1}$ during stable atmospheric conditions over non-vegetated areas (i.e., over bare ground or when snow entirely buries the vegetation), as follows in [energBalNoVegSolve.f90](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/energBalNoVegSolve.f90):

$$
Q_H = \left(\rho_{\text{air}}  c_P  C_H  U + E_0\right) \left(T_s - \theta_a\right),
$$

where $\rho_{\text{air}}$ is the density of air (kg m $^{-3}$), $c_P$ is the specific heat capacity of air (set to $1.00464 \dot 10^3$ J kg $^{-1}$ K $^{-1}$ in CLASSIC), $C_H$ is the surface drag coefficient (unitless), $U$ the wind speed at reference height (m s $^{-1}$), $T_s$ the surface temperature (K), and $\theta_a$ the potential air temperature at the reference height (K).  $E_0$ is set to 2 W m $^{-2}$ K $^{-1}$ when $T_s < \theta_a$ (i.e., for atmospheric stable condition) and 0 W m $^{-2}$ K $^{-1}$ otherwise.

See Brown et al. ([2006](https://doi.org/10.3137/ao.440302)) Fig. 11 for more details.

- Set $E_0 = 2$ W m $^{-2}$ K $^{-1}$ if there is snow: https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blame/SnowC2-1D/src/energBalNoVegSolve.f90#L269
- Add it in the sensible heat flux: e.g., https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/energBalNoVegSolve.f90#L594 (already coded by Ross Brown; see paper above)


---

### 2. Arctic Snowpack Adaptations

#### 2.1. Blowing snow sublimation losses 

The blowing snow sublimation loss parameterization of Gordon et al. ([2006](https://doi.org/10.3137/ao.440303)) is implemented:

$$
Q_s = 0.0018
\left(\frac{T_0}{T_{\text{air}}}\right)^4 
U_t \rho_a q_{si} (1-RH_i)
\left(\frac{U_{10}}{U_t}\right)^{3.6}, 
\quad \text{for } U_{10} > U_t  \text{ and }  T_{\text{air}} < T_0,
$$

with:

$$
U_t = 6.98 + 0.0033 \left(T_{\text{air}} - 245.88\right)^2,
$$

where $U_t$ [m s $^{-1}$] is the threshold wind speed at 10 m for the initiation of blowing snow, $T_{\text{air}}$ is the near-surface air temperature [K], $T_0$ is the freezing point of water (defined as 273.16 °K in CLASSIC), $\rho_a$ is the density of air [kg m $^{-3}$], $q_{si}$ is the saturation specific humidity of ice at reference height [kg kg $^{-1}$], $RH_i$ is the relative humidity with respect to ice [fraction], $U_{10}$ is the wind speed at 10 m above the snow surface [m s $^{-1}$]. This equation is only applied over bare ground subareas covered with snow (i.e., over bare ground and when snow buries the vegetation). The blowing snow sublimation losses from the intercepted snow on the canopy are not considered. 

- Implementation in a new file [snowWindSublimation.f90](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/snowWindSublimation.f90)
- 10 m wind speed adjustment in [energyBudgetDriver.f90](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/energyBudgetDriver.f90#L1153) (L1153)
- Blowing snow sublimation loss rate added to the evaporation flux in [energBalNoVegSolve.f90](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/energBalNoVegSolve.f90): 
  - https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/energBalNoVegSolve.f90#L600
  - https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/energBalNoVegSolve.f90#L902



---


#### 2.2. Snow compaction scheme  

##### Fresh snow density (not modified) in [atmosphericVarsCalc.f90](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/atmosphericVarsCalc.f90#L135)

##### Maximum snow density

The snow density $\rho_s$  in CLASSIC increases towards a maximum snow density $\rho_{\max }$ in an exponential way as follows:

$$
\rho_s(t+1) = \left[\rho_s(t) - \rho_{\max}\right] 
\exp  \left(-\frac{0.01 \Delta t}{3600}\right) + \rho_{\max},
$$

with:

$$
\begin{aligned}
\rho_{\max} &= 450 - \frac{204.7}{d_s}\left[1.0 - \exp  \left(-\frac{d_s}{0.673}\right)\right]
\quad \text{for } T_s < 0^{\circ}\mathrm{C}, \\
\rho_{\max} &= 700 - \frac{204.7}{d_s}\left[1.0 - \exp  \left(-\frac{d_s}{0.673}\right)\right]
\quad \text{for } T_s = 0^{\circ}\mathrm{C}.
\end{aligned}
$$


where $d_s$ is the snow depth [m], $T_s$ the snow temperature [°C], and $\Delta t$ the model time step [s]. The maximum snow density formulations are based on Tabler et al. ([1990](https://westernsnowconference.org/node/620)), with additional temperature threshold and enhanced settling rates determined empirically (Brown, unpublished manuscript, [2001](https://www.dropbox.com/scl/fo/iu8l96afrllug7qs24qwu/ALwo3FkohrOK6Pl52W8ULDg?rlkey=bh5od4983961yef5x99zzpugk&st=j2xsh20z&dl=0)).

The maximum snow density for dry snow was modified as follows depending on wind speed:
$$
\rho_{\max} =
430 - \frac{204.7}{d_s}
\left[1.0 - \exp  \left(-\frac{d_s}{0.673}\right)\right]
\quad \text{for } T_s < 0^{\circ}\mathrm{C}
\text{ and } U < 2.5~\mathrm{m s^{-1}}
\quad ( = \rho_{\text{no wind}} )
$$

$$
\begin{aligned}
\rho_{\max} &=
\rho_{\text{no wind}}
+ (\rho_{\text{wind}} - \rho_{\text{no wind}})
\exp  \left[-\frac{(d_s - d_0)^2}{2\sigma^2}\right]
\quad \text{for } T_s < 0^{\circ}\mathrm{C},\ U \ge 2.5~\mathrm{m s^{-1}} \\
&=
430 - \frac{204.7}{d_s}
\left[1.0 - \exp  \left(-\frac{d_s}{0.673}\right)\right]
\left\{
1 -
\left[1 - \exp  \left(-\frac{U}{U_0}\right)\right]
\exp  \left[-\frac{(d_s - d_0)^2}{2\sigma^2}\right]
\right\}
\end{aligned}
$$

with

$$
\rho_{\text{wind}} =
430 - \frac{204.7}{d_s}
\left[1.0 - \exp  \left(-\frac{d_s}{0.673}\right)\right]
\exp  \left(-\frac{U}{U_0}\right),
$$

where $U$ is the wind speed at 2 m above the snowpack [m s $^{-1}$], $U_0=2.5$ m s $^{-1}$ or $3.5$ m s $^{-1}$ (TVC adjusted), $d_0=0.8$ m, and $\sigma=1.0$ m. The constant $450$ is decreased to $430$ in order to better fit with the observations at all sites, considering the induced increased wind compaction during strong wind events. The wind induced compaction is used only over the bare ground subfraction (i.e., over non-vegetated subareas and/or when the snow buries the vegetation) as we consider that there is no further snow compaction within the canopy. Over the vegetated subarea $\rho_{\text {no wind }}$ is maintained regardless of the wind speed. We adjust the wind speed provided at specific measurement heights to 2 meters in the model using a logarithmic wind profile. For wet snow, the Brown et al. ([2006](https://doi.org/10.3137/ao.440302))'s equation is kept unchanged.

- Implementation of the new maximum snow density: https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/snowProcesses.f90#L167
- 2 m wind speed adjustment: https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/waterBudgetDriver.f90#L653

---


#### 2.3. Snow thermal conductivity

The current Sturm et al. ([1997](https://www.cambridge.org/core/journals/journal-of-glaciology/article/thermal-conductivity-of-seasonal-snow/7888DCB1F06AFFC755B6D4D567833925)) snow thermal conductivity parameterization:  

$$
\lambda_s =
\begin{cases}
0.234 \times 10^{-3} \rho_s + 0.023,
& \rho_s < 156~\mathrm{kg m^{-3}}, \\[0.5ex]
3.233 \times 10^{-6} \rho_s^2
- 1.01 \times 10^{-3} \rho_s + 0.138,
& \rho_s \ge 156~\mathrm{kg m^{-3}}.
\end{cases}
$$

was replaced with the Calonne et al. ([2011](https://onlinelibrary.wiley.com/doi/abs/10.1029/2011GL049234))'s one:

$$
\lambda_s =
2.5 \times 10^{-6} \rho_s^2
- 1.23 \times 10^{-4} \rho_s
+ 0.024,
$$

with $\lambda_s$ the snow thermal conductivity [W m $^{-1}$ K $^{-1}$], and $\rho_s$ the snow density [kg m $^{-3}$].

- Implementation of the new snow thermal conductivity: https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/src/energyBudgetPrep.f90#L666





## Model Testing and Evaluation

The modified model was tested and evaluated in **1D site simulations** at a set of mid-latitude, alpine, maritime, taiga, and Arctic sites.  

| Site                   | Short | Latitude  | Longitude  | Elevation | Period used | Type     | Reference                                 |
|------------------------|-------|-----------|------------|-----------|-------------|----------|-------------------------------------------|
| Col de Porte, France   | cdp   | 45.30° N  | 5.77° E    | 1325 m    | 1994–2014   | Alpine   | Morin et al. ([2012](https://essd.copernicus.org/articles/4/13/2012/)), Lejeune et al. ([2019](https://essd.copernicus.org/articles/11/71/2019/)) |
| Reynolds Mountain East, USA | rme | 43.06° N  | 116.75° W | 2060 m    | 1988–2008   | Alpine   | Reba et al. ([2011](https://onlinelibrary.wiley.com/doi/abs/10.1029/2010WR010030)) |
| Senator Beck, USA      | snb   | 37.91° N  | 107.73° W  | 3714 m    | 2005–2015   | Alpine   | Landry et al. ([2014](https://onlinelibrary.wiley.com/doi/abs/10.1002/2013WR013711)) |
| Swamp Angel, USA       | swa   | 37.91° N  | 107.71° W  | 3371 m    | 2005–2015   | Alpine   | Landry et al. ([2014](https://onlinelibrary.wiley.com/doi/abs/10.1002/2013WR013711)) |
| Sapporo, Japan         | sap   | 43.08° N  | 141.34° E  | 15 m      | 2005–2015   | Maritime | Niwano et al. ([2012](https://onlinelibrary.wiley.com/doi/abs/10.1029/2011JF002239)) |
| Sodankylä, Finland     | sod   | 67.37° N  | 26.63° E   | 179 m     | 2007–2014   | Taiga    | Essery et al. ([2016](https://gi.copernicus.org/articles/5/219/2016/)) |
| Weissfluhjoch, Switzerland | wfj | 46.83° N  | 9.81° E    | 2540 m    | 1996–2016   | Alpine   | Wever et al. ([2015](https://tc.copernicus.org/articles/9/2271/2015/)), WSL ([2017](https://www.envidat.ch/dataset/snowmip)) |
| Bylot Island, Canada   | byl   | 73.15° N  | 80.00° W   | 25 m      | 2014–2019   | Arctic   | Domine et al. ([2021](https://essd.copernicus.org/articles/13/4331/2021/)) |
| Umiujaq TUNDRA, Canada | umt   | 56.56° N  | 76.48° W   | 132 m     | 2016–2021   | Arctic   | Domine et al. ([2024](https://essd.copernicus.org/articles/16/1523/2024/)) |
| Trail Valley Creek, Canada | tvc | 56.55° N  | 76.47° W   | 82 m      | 2017–2019   | Arctic   | Dutch et al. ([2022](https://tc.copernicus.org/articles/16/4201/2022/)), Boike et al. ([2023](https://doi.pangaea.de/10.1594/PANGAEA.962726)) |



## Experiments

The following experiments were run as part of the study:

| Experiment | Description | Binary file |
|------------|-------------|-------------|
| DEF        | Default model version | `DEF` |
| TCZERO     | DEF with the thermal conductivity bug correction | `BUG_CORRECT` |
| TSNBOT     | TCZERO with improved bottom snow temperature | `BUG_CORRECT_TSNBT_OP1` |
| EZERO (PHYS) | TSNBOT with windless exchange coefficient (all physics improvements) | `BUG_CORRECT_TSNBT_OP1_EZERO` |
| SUBLI      | PHYS with blowing snow sublimation losses | `PHYS_ALL_SUBLI_v2` |
| COMPAC     | SUBLI with improved snow compaction scheme | `PHYS_ALL_SUBLI_v2_COMPAC_v1` |
| CL11 (ALL) | COMPAC with Calonne et al. ([2011](https://onlinelibrary.wiley.com/doi/abs/10.1029/2011GL049234)) snow conductivity parameterization (all developments) | `PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne` |
| ALL_TVC    | ALL with adjusted max snow density parameters including TVC | `PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne` |



## Job Options Files

Each site and experiment has an associated job options file, which specifies the forcing, initialization, and output directories.  



### SnowMIP Sites
Base path:  
`https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/<site>/job_options_file_run_<experiment>.txt`

Experiments included per site:  
- DEF  
- BUG_CORRECT_TSNBT_OP1_EZERO  
- PHYS_ALL_SUBLI_v2  
- PHYS_ALL_SUBLI_v2_COMPAC_v1  
- PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne  
- PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne  

**Col de Porte (cdp)**  
- [DEF](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/cdp/job_options_file_run_DEF.txt)  
- [BUG_CORRECT_TSNBT_OP1_EZERO](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/cdp/job_options_file_run_BUG_CORRECT_TSNBT_OP1_EZERO.txt)  
- [PHYS_ALL_SUBLI_v2](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/cdp/job_options_file_run_PHYS_ALL_SUBLI_v2.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/cdp/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/cdp/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/cdp/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne.txt)  

**Reynolds Mountain East (rme)**  
- [DEF](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/rme/job_options_file_run_DEF.txt)  
- [BUG_CORRECT_TSNBT_OP1_EZERO](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/rme/job_options_file_run_BUG_CORRECT_TSNBT_OP1_EZERO.txt)  
- [PHYS_ALL_SUBLI_v2](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/rme/job_options_file_run_PHYS_ALL_SUBLI_v2.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/rme/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/rme/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/rme/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne.txt)  

**Senator Beck (snb)**  
- [DEF](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/snb/job_options_file_run_DEF.txt)  
- [BUG_CORRECT_TSNBT_OP1_EZERO](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/snb/job_options_file_run_BUG_CORRECT_TSNBT_OP1_EZERO.txt)  
- [PHYS_ALL_SUBLI_v2](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/snb/job_options_file_run_PHYS_ALL_SUBLI_v2.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/snb/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/snb/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/snb/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne.txt)  

**Swamp Angel (swa)**  
- [DEF](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/swa/job_options_file_run_DEF.txt)  
- [BUG_CORRECT_TSNBT_OP1_EZERO](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/swa/job_options_file_run_BUG_CORRECT_TSNBT_OP1_EZERO.txt)  
- [PHYS_ALL_SUBLI_v2](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/swa/job_options_file_run_PHYS_ALL_SUBLI_v2.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/swa/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/swa/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/swa/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne.txt)  

**Sapporo (sap)**  
- [DEF](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sap/job_options_file_run_DEF.txt)  
- [BUG_CORRECT_TSNBT_OP1_EZERO](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sap/job_options_file_run_BUG_CORRECT_TSNBT_OP1_EZERO.txt)  
- [PHYS_ALL_SUBLI_v2](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sap/job_options_file_run_PHYS_ALL_SUBLI_v2.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sap/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sap/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sap/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne.txt)  

**Sodankylä (sod)**  
- [DEF](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sod/job_options_file_run_DEF.txt)  
- [BUG_CORRECT_TSNBT_OP1_EZERO](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sod/job_options_file_run_BUG_CORRECT_TSNBT_OP1_EZERO.txt)  
- [PHYS_ALL_SUBLI_v2](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sod/job_options_file_run_PHYS_ALL_SUBLI_v2.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sod/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sod/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/sod/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne.txt)  

**Weissfluhjoch (wfj)**  
- [DEF](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/wfj/job_options_file_run_DEF.txt)  
- [BUG_CORRECT_TSNBT_OP1_EZERO](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/wfj/job_options_file_run_BUG_CORRECT_TSNBT_OP1_EZERO.txt)  
- [PHYS_ALL_SUBLI_v2](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/wfj/job_options_file_run_PHYS_ALL_SUBLI_v2.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/wfj/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/wfj/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowMIP/wfj/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne.txt)  




### Arctic Sites
Base path:  
`https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/<site>/job_options_file_run_<experiment>.txt`

**Bylot Island (byl)**  
- [peat_DEF](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/byl/job_options_file_run_peat_DEF.txt)  
- [peat_BUG_CORRECT_TSNBT_OP1_EZERO](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/byl/job_options_file_run_peat_BUG_CORRECT_TSNBT_OP1_EZERO.txt)  
- [peat_PHYS_ALL_SUBLI_v2](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/byl/job_options_file_run_peat_PHYS_ALL_SUBLI_v2.txt)  
- [peat_PHYS_ALL_SUBLI_v2_COMPAC_v1](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/byl/job_options_file_run_peat_PHYS_ALL_SUBLI_v2_COMPAC_v1.txt)  
- [peat_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/byl/job_options_file_run_peat_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne.txt)  
- [peat_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/byl/job_options_file_run_peat_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne.txt)  

**Umiujaq (umt)**  
- [DEF_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/umt/job_options_file_run_DEF_correct_SH.txt)  
- [BUG_CORRECT_TSNBT_OP1_EZERO_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/umt/job_options_file_run_BUG_CORRECT_TSNBT_OP1_EZERO_correct_SH.txt)  
- [PHYS_ALL_SUBLI_v2_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/umt/job_options_file_run_PHYS_ALL_SUBLI_v2_correct_SH.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/umt/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1_correct_SH.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/umt/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne_correct_SH.txt)  
- [PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/umt/job_options_file_run_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne_correct_SH.txt)  

**Trail Valley Creek (tvc)**  
- [1peat_2xSnowf_DEF_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/tvc/job_options_file_run_1peat_2xSnowf_DEF_correct_SH.txt)  
- [1peat_2xSnowf_BUG_CORRECT_TSNBT_OP1_EZERO_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/tvc/job_options_file_run_1peat_2xSnowf_BUG_CORRECT_TSNBT_OP1_EZERO_correct_SH.txt)  
- [1peat_2xSnowf_PHYS_ALL_SUBLI_v2_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/tvc/job_options_file_run_1peat_2xSnowf_PHYS_ALL_SUBLI_v2_correct_SH.txt)  
- [1peat_2xSnowf_PHYS_ALL_SUBLI_v2_COMPAC_v1_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/tvc/job_options_file_run_1peat_2xSnowf_PHYS_ALL_SUBLI_v2_COMPAC_v1_correct_SH.txt)  
- [1peat_2xSnowf_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/tvc/job_options_file_run_1peat_2xSnowf_PHYS_ALL_SUBLI_v2_COMPAC_v1_calonne_correct_SH.txt)  
- [1peat_2xSnowf_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne_correct_SH](https://github.com/mickaellalande/CLASSIC-SnowC2-1D/blob/SnowC2-1D/inputFiles/SnowArctic/tvc/job_options_file_run_1peat_2xSnowf_PHYS_ALL_SUBLI_v2_COMPAC_2.5_LIM_3.5_calonne_correct_SH.txt)  



Note: Additional testing experiments and site simulations are also included in this repository, beyond those described in the published paper.




## Citation

If you use this version of the model in your work, please cite the corresponding paper (citation details will be updated once the GMD article is published).



## Disclaimer

This repository is a **research copy** intended solely to accompany the publication listed above.  
It should **not** be considered the official or maintained version of CLASSIC.
