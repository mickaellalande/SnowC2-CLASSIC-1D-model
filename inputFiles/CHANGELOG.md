# Changelog

The numbering here should use the same numbering as the main code for release versions but minor increments need not match up.

## 1.0.10

Add new N cycle variables to init-files

## 1.0.9 

### Correct vegetation composition and soil textures at some sites

- RU-SkP: This Larch forest site was previously set up as 100% NdlDcdTr, but has a lower stand density and a shrub understory. Now set up as 72% NdlDcdTr and 28% BdlEvgSh.
- ES-LgS: The site was previously set up as 45% GrassC3, but it consists of both grasses and shrubs. Now set up as 22% GrassC3 and 23% BdlEvgSh.
- IT-Noe: The soil permeable depth was increased from 0.4 to 1.0 m, DRN set to 0 and the soil texture was modified.
- US-Whs: Modify the soil texture and set DRN = 0.

## 1.0.8
### Add in the Samoylov (Lena River Delta) (RU-Sam) site (Jan 13 2023)

Adding in this site setup by Rose Lefebvre

## 1.0.7 

### Add in the Trail Valley Creek dwarf-shrub tundra site (CA-TVC) (Nov 21 2022)
## 1.0.6

### Addition of new variables (Jan 5 2022)

- Added new variables required to fix the restart bug.
- Added competition variables to init-files that did not include them yet.
## 1.0.5

### Add in 8 peatland sites originally published in Wu et al. GMD (2016)

- Sites include CA-Mer, CA-WP1

## 1.0.4 

### Add in 7 shrub sites from the FLUXNET2015 dataset (DK-ZaH, ES-Amo, ES-LJu, IT-Noe, US-SRC, US-Sta, US-Whs) and the Daring Lake dwarf-shrub tundra site (CA-DL1)

### Correct setup for the Siberian deciduous needleleaf forest site RU-SkP 

- The site is a Larch forest, but was previously set up as 100% DCD-CLD BdlTr. Now it is changed to 100% DCD NdlTr.

## 1.0.3

### Add in the temperate deciduous broadleaf sites from Asaadi et al. (2018): US-MMS, US-Ha1, US-UMB

- To check on the non-structural carbohydrates parameterization we need to have the US-MMS, US-Ha1, US-UMB sites integrated into the benchmarking. 
- Note that the sites are not setup exactly the same if better information was found or errors were corrected. See each sites yaml file for any relevant info.

## 1.0.2

### Expand PFT arrays from 9 PFTs to 12 (Feb 9 2021)

- The basic array of 9 PFTs is expanded to include two shrubs and sedges.
- The present benchmarking suite does not make use of these open positions but it does allow the use of a standard parameter file that includes them.
- Edited the yaml files to allow information about leap years in the meteorology and the PFTs that are present at the site.

## 1.0.1

### Addition of new variables (Feb 9 2021)

- Added non-structural and structural carbon pools to leaves, roots, and stems.
- Added peat soil C pool
- Convert litter and soil C to be tracked per layer

## 1.0

### Initial release

- Initial release associated with CLASSIC 1.0. 
- Available from [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3525336.svg)](https://doi.org/10.5281/zenodo.3525336)
