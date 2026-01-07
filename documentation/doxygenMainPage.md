CLASSIC main page {#mainpage}
============

# The Canadian Land Surface Scheme Including biogeochemical Cycles (CLASSIC) {#main}

The Canadian Land Surface Scheme Including biogeochemical Cycles (CLASSIC) simulates the exchanges of energy, water, carbon, and momentum at the earth's surface. CLASSIC is formed by the coupling of the Canadian Land Surface Scheme (CLASS) and the Canadian Terrestrial Ecosystem Model (CTEM). CLASS handles the model physics including fluxes of energy, water and momentum. CTEM simulates biogeochemical cycles including fluxes of carbon and nitrogen.

# [Link to the CLASSIC webpage](https://cccma.gitlab.io/classic_pages)

\image html CLASSIC_run_schematic.png "CLASSIC code structure" width=850px

1. @subpage overviewCLASSIC
   - @subpage modelVers
   - @subpage PFTsCLASSIC
   - @subpage overviewCLASS "Historical overview of the Canadian Land Surface Scheme (CLASS)"
     - @subpage devHistory
   - @subpage overviewCTEM "Historical overview of the Canadian Terrestrial Ecosystem Model (CTEM)"
2. @subpage tiles "Tiling for sub-grid heterogeneity"
  - @subpage compvsmosaic
  - @subpage dynTile
3. @subpage basicInputs "The five basic types of data that are required to run CLASSIC"
  - @subpage modelParams
  - @subpage forcingData
  - @subpage vegetationData
    - @subpage vegCLASSonly
    - @subpage vegCTEMtoo
  - @subpage soilData
  - @subpage initProgVar
    - @subpage initPhysProgVar
    - @subpage initCTEMProgVar
4. @subpage CTEMaddInputs
  - Greenhouse gases
    - @subpage initCO2
    - @subpage initCH4
  - Disturbance (fire) inputs
    - @subpage initLightFire
    - @subpage initPopd
  - @subpage initClimComp "Competition for space between PFTs"
  - Prognostic simulation of methane emissions
    - @subpage initWetSlope "Dynamically-determined wetlands"
    - @subpage initWetArea
  - @subpage initPeat "Peatlands"
  - @subpage inputLUC
  - @subpage tracers "Carbon tracers"
  - @subpage nitroData "Nitrogen cycle"
5. @subpage makeInputFiles
  - @subpage makeMet "Meteorological inputs"
  - @subpage makeInit "The model initialization and restart files"
  - @subpage ghgfiles "Greenhouse gas inputs"
  - @subpage makeOther "Other inputs"
6. @subpage runPrep "Preparing a CLASSIC run"
  - @subpage Environ
      - @subpage Containers
  - @subpage compilingMod
  - @subpage setupJobOpts
  - @subpage xmlSystem "Configuring the model outputs"
7. @subpage runCLASSIC "Running CLASSIC"
  - @subpage runStandAloneMode "Running CLASSIC for a point location"
  - @subpage runGrid "Running CLASSIC over a grid (regional,global)"
8. @subpage benchAmber
9. Tools available   
  - @subpage makeGHGfiles
  - @subpage xmlSystem
  - @subpage asciiMet
  - @subpage initTool
  - @subpage modifyRS
  - @subpage convertToCSV
  - @subpage regTest
10. @subpage howDoI
11. @subpage legacyFileNames
12. [CLASSIC Code Conventions](http://cccma.gitlab.io/classic_pages/info/conventions/)
