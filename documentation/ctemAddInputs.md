# Additional inputs depending on model configuration  {#CTEMaddInputs}

CLASSIC can be run in several different configurations, some of which require additional model inputs including:

1. Greenhouse gases
  1. @ref initCO2
  2. @ref initCH4
2. Disturbance (fire) inputs
  1. @ref initLightFire
  2. @ref initPopd
3. [Competition for space between PFTs](@ref initClimComp)
4. Prognostic simulation of methane emissions
  1. [Dynamically-determined wetlands](@ref initWetSlope)
  2.  @ref initWetArea
5. [Peatlands](@ref initPeat)
6. @ref inputLUC
7. [Tracers](@ref tracers)
8. [Nitrogen cycle](@ref Ncycle)

---

# Atmospheric carbon dioxide concentration {#initCO2}

**Annual** atmospheric carbon dioxide concentrations are needed for CLASS+CTEM runs. The annual values are read in from a netcdf file. If you use the NCO tool ncdump (included in CLASSIC Singularity container), and do `ncdump -hs` on a properly formatted file you should yield something similar to below (only relevant sections shown here for a file with 321 years of data). Note the time units. 

        [scrd530@ppp5login-002 CO2]$ ncdump -h TRENDY_CO2_1700-2020_GCP2021.nc
		netcdf TRENDY_CO2_1700-2020_GCP2021 {
		dimensions:
				time = 321 ;
		variables:
				double time(time) ;
						time:standard_name = "time" ;
						time:units = "day as YYYYMMDD.FFFF" ;
						time:calendar = "proleptic_gregorian" ;
				double CO2(time) ;
						CO2:_FillValue = 1.e+38 ;
						CO2:title = "Atmospheric Carbon Dioxide Concentration" ;
						CO2:units = "ppmv" ;

 The variable name is not important as long as it is the only variable in the file besides time. Units expected are ppmv.

# Atmospheric methane concentration {#initCH4}

**Annual** atmospheric methane concentrations are needed for CLASSIC runs with biogeochemistry on (CTEM on). The annual values are read in from a netcdf file similar to CO2. The file format is the same as CO2. The variable name is not important as long as it is the only variable in the file besides time. Units expected are ppmv.

# Lightning frequency for fire ignition {#initLightFire}

**Daily cloud-to-ground** lightning frequency is used by the disturbance subroutine for fire. The code at present is set to use daily lightning frequency. If you have mean monthly values you can interpolate them to daily values (`cdo inttime` is useful here). An `ncdump -hs` of a properly formatted file is below. Note the file is chunked for a T63 grid (128 x 64), other grids may require different chunk sizes for optimal performance. The variable name is not important as long as it is the only variable in the file besides lat ,lon, and time. Note the units of the lght_lisotd and time variables.

        netcdf lisotd_1995_2014_climtlgl_lghtng_as_ts_1700_2050_chunked {
        dimensions:
        	lat = 64 ;
        	time = UNLIMITED ; // (128535 currently)
        	lon = 128 ;
        variables:
        	double lat(lat) ;
        		lat:standard_name = "latitude" ;
        		lat:long_name = "latitude" ;
        		lat:units = "degrees_north" ;
        		lat:axis = "Y" ;
        		lat:_Storage = "contiguous" ;
        		lat:_Endianness = "little" ;
        	float lght_lisotd(time, lat, lon) ;
        		lght_lisotd:long_name = "Combined C2G (see Other_info in global attributes) Flash Rate Annual Climatology (1995-2014)" ;
        		lght_lisotd:units = "**strikes km-2 yr-1**" ;
        		lght_lisotd:grid_type = "gaussian" ;
        		lght_lisotd:_FillValue = -1.e+38f ;
        		lght_lisotd:missing_value = -1.e+38f ;
        		lght_lisotd:_Storage = "chunked" ;
        		lght_lisotd:_ChunkSizes = 128535, 8, 16 ;
        		lght_lisotd:_Endianness = "little" ;
        	double lon(lon) ;
        		lon:standard_name = "longitude" ;
        		lon:long_name = "longitude" ;
        		lon:units = "degrees_east" ;
        		lon:axis = "X" ;
        		lon:_Storage = "contiguous" ;
        		lon:_Endianness = "little" ;
        	double time(time) ;
        		time:standard_name = "time" ;
        		time:units = "**day as %Y%m%d.%f**" ;
        		time:calendar = "standard" ;
        		time:_Storage = "chunked" ;
        		time:_ChunkSizes = 128535 ;
        		time:_Endianness = "little" ;


# Population density for fire ignition/suppresion {#initPopd}

Fire uses a time series of **annually** varying population density for fire suppression and ignition. An `ncdump -hs` of a properly formatted file is below. Note the file is chunked for a T63 grid (128 x 64), other grids may require different chunk sizes for optimal performance. The variable name is not important as long as it is the only variable in the file besides lat ,lon, and time.

          netcdf POPD_annual_1700_2017_T63_chunked {
          dimensions:
          	lat = 64 ;
          	lon = 128 ;
          	time = UNLIMITED ; // (318 currently)
          variables:
          	double lat(lat) ;
          		lat:standard_name = "latitude" ;
          		lat:long_name = "latitude" ;
          		lat:units = "degrees_north" ;
          		lat:axis = "Y" ;
          		lat:_Storage = "contiguous" ;
          		lat:_Endianness = "little" ;
          	double lon(lon) ;
          		lon:standard_name = "longitude" ;
          		lon:long_name = "longitude" ;
          		lon:units = "degrees_east" ;
          		lon:axis = "X" ;
          		lon:_Storage = "contiguous" ;
          		lon:_Endianness = "little" ;
          	float popd(time, lat, lon) ;
          		popd:grid_type = "gaussian" ;
              popd:units = **"Number of people / km2"** ;
          		popd:_FillValue = -9999.f ;
          		popd:missing_value = -9999.f ;
          		popd:_Storage = "chunked" ;
          		popd:_ChunkSizes = 318, 8, 16 ;
          		popd:_Endianness = "little" ;
          	double time(time) ;
          		time:standard_name = "time" ;
          		time:units = "**day as %Y%m%d.%f**" ;
          		time:calendar = "proleptic_gregorian" ;
          		time:_Storage = "chunked" ;
          		time:_ChunkSizes = 318 ;
          		time:_Endianness = "little" ;

# Climatic variables for PFT competition simulations {#initClimComp}

The PFT competition scheme ( @ref competitionMod.f90 ) uses bioclimatic variables to determine if a PFT can establish and attempt to colonize a grid cell (see competitionMod.f90:bioclim and competitionMod.f90:existence). These climatic variables are written to the model restart file for use to initialize the next run. They are also either read in from the initialization file or are determined during a model run (see model switch *inibioclim* in the job options file). They can be set to values of 0 in the initial conditions file and then allow *inibioclim* to properly initialize them.

- *anndefct* Annual water deficit, i.e. daily values of potential evaporation that exceed precipitation accumulated over a year [mm]
- *annsrpls* Annual water deficit, i.e. daily values of precipitation that exceed potential evaporation accumulated over a year [mm]
- *annpcp* Annual precipitation [mm]
- *aridity* Aridity index, ratio of potential evaporation to precipitation [ ]
- *defctmon* Number of months in a year with water deficit, i.e. precipitation less than potential evaporation [months]
- *srplsmon* Number of months in a year with surplus water, i.e. precipitation more than potential evaporation [months]
- *dry_season_length* Length of consecutive dry season in months, where a dry month is defined as the month in which potential evaporation exceeds precipitation [months]
- *gdd5* Growing degree days above 5 C [days]
- *tcoldm* Temperature of the coldest month [deg C]
- *twarmm* Temperature of the warmest month [deg C]

# Orographic information for dynamic wetland scheme {#initWetSlope}

Eight slope based fractions are read in from the model initialization file for calculating dynamic wetland fractions. As the soil moisture in a grid cell increases above specified thresholds then the really flat portions of the grid cell are assumed to gradually turn into wetlands. The eight slope based fractions correspond to the fraction of the grid cell that have slope less than 0.025%, 0.05%, 0.1%, 0.15%, 0.20%, 0.25%, 0.3% and 0.35%. The relevant variables in the initialization file are shown below. In the file, *slope* is a dimension and *slopefrac* is a variable.

        double slope(slope) ;
          slope:long_name = "wetland slope fractions for 0.025, 0.05, 0.1, 0.15, 0.20, 0.25, 0.3 and 0.35 percent slope threshold" ;
          ...
        float slopefrac(slope, tile, lat, lon) ;
          slopefrac:_FillValue = -999.f ;
          slopefrac:units = "-" ;
          slopefrac:long_name = "Slope-based fraction for dynamic wetlands" ;

This formulation is sensible regionally but shouldn't be used at the site-level.
# Prescribed wetland area {#initWetArea}

**Daily** values of wetland fraction are used for modelling methane emissions from wetlands.

          netcdf gcp-ch4_wetlands_1838-2017_t63_final_daily {
          dimensions:
                  time = UNLIMITED ; // (65700 currently)
                  lat = 64 ;
                  lon = 128 ;
          variables:
                  double Fw(time, lat, lon) ;
                          Fw:long_name = "Fraction inundated" ;
                          Fw:units = "**fraction**" ;
                          Fw:grid_type = "gaussian" ;
                          Fw:_FillValue = -9999. ;
                          Fw:missing_value = -9999. ;
                          Fw:_Storage = "chunked" ;
                          Fw:_ChunkSizes = 65700, 8, 16 ;
                  double lat(lat) ;
                          lat:standard_name = "latitude" ;
                          lat:long_name = "latitude" ;
                          lat:units = "degrees_north" ;
                          lat:axis = "Y" ;
                          lat:_Storage = "contiguous" ;
                  double lon(lon) ;
                          lon:standard_name = "longitude" ;
                          lon:long_name = "longitude" ;
                          lon:units = "degrees_east" ;
                          lon:axis = "X" ;
                          lon:_Storage = "contiguous" ;
                  double time(time) ;
                          time:standard_name = "time" ;
                          time:long_name = "time" ;
                          time:units = "**day as %Y%m%d.%f**" ;
                          time:calendar = "365_day" ;
                          time:_Storage = "chunked" ;
                          time:_ChunkSizes = 65700 ;



# Peatland variables {#initPeat}

Peatlands are simulated following the parameterization of Wu et al. (2016) \cite Wu2016-zt. Please read about them [here](@ref peatlandC). The peatland areas are specified by the *ipeatland* flag in the initialization file:

        float ipeatland(tile, lat, lon) ;
          ipeatland:_FillValue = -999.f ;
          ipeatland:units = "-" ;
          ipeatland:long_name = "Peatland flag: 0 = not a peatland, 1= bog, 2 = fen" ;

There are several prognostic variables that are associated with the peatland areas of the gridcell. These may be initialized to zero prior to a spinup.

        float Cmossmas(tile, lat, lon) ;
          Cmossmas:_FillValue = -999.f ;
          Cmossmas:units = "kgC/m2" ;
          Cmossmas:long_name = "C in moss biomass" ;

        float litrmsmoss(tile, lat, lon) ;
          litrmsmoss:_FillValue = -999.f ;
          litrmsmoss:units = "kgC/m2" ;
          litrmsmoss:long_name = "Moss litter mass" ;

        float dmoss(tile, lat, lon) ;
          dmoss:_FillValue = -999.f ;
          dmoss:units = "m" ;
          dmoss:long_name = "Depth of living moss" ;

        float peatSoilC(tile, lat, lon) ;
		  peatSoilC:_FillValue = -999.f ;
		  peatSoilC:long_name = "Peat soil C for peatland tiles" ;
		  peatSoilC:units = "kg C/m^2" ;		

If you are running a single site peatland then the peatland tile is the whole grid cell. If you are running large regions you may wish to have the peatlands as a separate tile. This is done by having an *nmtest* > 1 and setting up each tile appropriately.

          int nmtest(lat, lon) ;
            nmtest:_FillValue = -999 ;
            nmtest:long_name = "Number of tiles in each grid cell" ;

# Land use change (LUC) {#inputLUC}

The LUC file contains an **annual** time series of fractional coverage of each of the CLASSIC PFTs. An `ncdump -hs` of a properly formatted file is below. Note the file is chunked for a T63 grid (128 x 64), other grids may require different chunk sizes for optimal performance. The variable name is not important as long as it is the only variable in the file besides lat ,lon, **lev** and time.

        netcdf GCP_2018_land_cover_CTEM_fractions_1700_2018_T63_chunked {
        dimensions:
        	time = UNLIMITED ; // (319 currently)
        	**lev = 9** ; ! This corresponds to number of biogeochemical (CTEM) PFTs.
        	lat = 64 ;
        	lon = 128 ;
        variables:
        	float frac(time, lev, lat, lon) ;
        		frac:grid_type = "gaussian" ;
            frac:units = **"fraction"** ;
        		frac:_FillValue = -9999.f ;
        		frac:missing_value = -9999.f ;
        		frac:_Storage = "chunked" ;
        		frac:_ChunkSizes = 319, 9, 8, 16 ;
        		frac:_Endianness = "little" ;
        	double lat(lat) ;
        		lat:standard_name = "latitude" ;
        		lat:long_name = "latitude" ;
        		lat:units = "degrees_north" ;
        		lat:axis = "Y" ;
        		lat:_Storage = "contiguous" ;
        		lat:_Endianness = "little" ;
        	double lev(lev) ;
        		lev:axis = "Z" ;
        		lev:_Storage = "contiguous" ;
        		lev:_Endianness = "little" ;
        	double lon(lon) ;
        		lon:standard_name = "longitude" ;
        		lon:long_name = "longitude" ;
        		lon:units = "degrees_east" ;
        		lon:axis = "X" ;
        		lon:_Storage = "contiguous" ;
        		lon:_Endianness = "little" ;
        	double time(time) ;
        		time:standard_name = "time" ;
        		time:units = **"day as %Y%m%d.%f"** ;
        		time:calendar = "proleptic_gregorian" ;
        		time:_Storage = "chunked" ;
        		time:_ChunkSizes = 319 ;
        		time:_Endianness = "little" ;

# Tracers {#tracers}

In development.

# Nitrogen cycle {#Ncycle}

The DEP file contains an **annual** time series of N deposition rate. An `ncdump -hs` of a properly formatted file is below. Note the file is chunked for a T63 grid (128 x 64), other grids may require different chunk sizes for optimal performance.

        netcdf t63_1700_2020_input4MIPs_hist_plus_ssp585_nitrogen_deposition {
        dimensions:
        	time = UNLIMITED ; // (321 currently)
        	lon = 128 ;
        	lat = 64 ;
        variables:
        double time(time) ;
                time:standard_name = "time" ;
                time:units = "day as %Y%m%d.%f" ;
                time:calendar = "proleptic_gregorian" ;
                time:axis = "T" ;
                time:_Storage = "chunked" ;
                time:_ChunkSizes = 512 ;
                time:_Endianness = "little" ;
                time:_NoFill = "true" ;
        float lon(lon) ;
                lon:standard_name = "longitude" ;
                lon:long_name = "longitude" ;
                lon:units = "degrees_east" ;
                lon:axis = "X" ;
                lon:_Storage = "contiguous" ;
                lon:_Endianness = "little" ;
                lon:_NoFill = "true" ;
        float lat(lat) ;
                lat:standard_name = "latitude" ;
                lat:long_name = "latitude" ;
                lat:units = "degrees_north" ;
                lat:axis = "Y" ;
                lat:_Storage = "contiguous" ;
                lat:_Endianness = "little" ;
                lat:_NoFill = "true" ;
        float Ndep(time, lat, lon) ;
                Ndep:standard_name = "WetDryNHxNOy_N_deposition" ;
                Ndep:long_name = "Wet and dry, NHx and NOy deposition " ;
                Ndep:units = "$ g N m^{-2}yr^{-1} $" ;
                Ndep:original_name = "Wet and dry NHx and NOy deposition" ;
                Ndep:cell_methods = "time:mean" ;
                Ndep:_Storage = "chunked" ;
                Ndep:_ChunkSizes = 1, 64, 128 ;
                Ndep:_Endianness = "little" ;
                Ndep:_NoFill = "true" ;

The FER file contains an **annual** time series of N fertiliser application rate. An `ncdump -hs` of a properly formatted file is below. Note the file is chunked for a T63 grid (128 x 64), other grids may require different chunk sizes for optimal performance.

        netcdf t63_NMIP_N_fertilizer_1700_2020 {
		dimensions:
			time = UNLIMITED ; // (321 currently)
			lon = 128 ;
			lat = 64 ;
		variables:
			double time(time) ;
					time:standard_name = "time" ;
					time:units = "day as %Y%m%d.%f" ;
					time:calendar = "proleptic_gregorian" ;
					time:axis = "T" ;
					time:_Storage = "chunked" ;
					time:_ChunkSizes = 512 ;
					time:_Endianness = "little" ;
					time:_NoFill = "true" ;
			float lon(lon) ;
					lon:standard_name = "longitude" ;
					lon:long_name = "longitude" ;
					lon:units = "degrees_east" ;
					lon:axis = "X" ;
					lon:_Storage = "contiguous" ;
					lon:_Endianness = "little" ;
					lon:_NoFill = "true" ;
			float lat(lat) ;
					lat:standard_name = "latitude" ;
					lat:long_name = "latitude" ;
					lat:units = "degrees_north" ;
					lat:axis = "Y" ;
					lat:_Storage = "contiguous" ;
					lat:_Endianness = "little" ;
					lat:_NoFill = "true" ;
			float nfer(time, lat, lon) ;
					nfer:standard_name = "fertilizer_application_rate" ;
					nfer:long_name = "fertilizer application rate" ;
					nfer:_FillValue = -9999.f ;
					nfer:missing_value = "-9999.0" ;
					nfer:units = "g N m-2 cropland yr-1" ;
					nfer:_Storage = "chunked" ;
					nfer:_ChunkSizes = 1, 64, 128 ;
					nfer:_Endianness = "little" ;
					nfer:_NoFill = "true" ;

There are several prognostic variables that are associated with the N cycle. ngleafmas_ns, ngleafmas_s, nbleafmas, 
nrootmass_ns, nrootmass_s, nstemmass_ns, nstemmass_s, soilnmas, nlitrmass, nh4_mass, and no3_mass may be initialised
following their C equivalents and corresponding C:N ratios.
bnf_nat and bnf_ant may be initialized to zero prior to a spinup.

        float ngleafmas_ns(tile, icc, lat, lon) ;
                ngleafmas_ns:_FillValue = NaNf ;
                ngleafmas_ns:long_name = "non-structural green leaf nitrogen mass" ;
                ngleafmas_ns:units = "g N/m2" ;
        float ngleafmas_s(tile, icc, lat, lon) ;
                ngleafmas_s:_FillValue = NaNf ;
                ngleafmas_s:long_name = "structural green leaf nitrogen mass" ;
                ngleafmas_s:units = "g N/m2" ;
        float nbleafmas(tile, icc, lat, lon) ;
                nbleafmas:_FillValue = NaNf ;
                nbleafmas:long_name = "Brown leaf nitrogen mass" ;
                nbleafmas:units = "g N/m2" ;
        float nrootmass_ns(tile, icc, lat, lon) ;
                nrootmass_ns:_FillValue = NaNf ;
                nrootmass_ns:long_name = "non-structural root nitrogen mass" ;
                nrootmass_ns:units = "g N/m2" ;
        float nrootmass_s(tile, icc, lat, lon) ;
                nrootmass_s:_FillValue = NaNf ;
                nrootmass_s:long_name = "structural root nitrogen mass" ;
                nrootmass_s:units = "g N/m2" ;
        float nstemmass_ns(tile, icc, lat, lon) ;
                nstemmass_ns:_FillValue = NaNf ;
                nstemmass_ns:long_name = "non-structural stem nitrogen mass" ;
                nstemmass_ns:units = "g N/m2" ;
        float nstemmass_s(tile, icc, lat, lon) ;
                nstemmass_s:_FillValue = NaNf ;
                nstemmass_s:long_name = "structural stem nitrogen mass" ;
                nstemmass_s:units = "g N/m2" ;
        float soilnmas(tile, iccp2, lat, lon) ;
                soilnmas:_FillValue = NaNf ;
                soilnmas:long_name = "Soil N mass per soil layer" ;
                soilnmas:units = "g N/m2" ;
        float nlitrmass(tile, iccp2, lat, lon) ;
                nlitrmass:_FillValue = NaNf ;
                nlitrmass:long_name = "Litter nitrogen mass per soil layer" ;
                nlitrmass:units = "g N/m2" ;
        float nh4_mass(tile, iccp2, lat, lon) ;
                nh4_mass:_FillValue = NaNf ;
                nh4_mass:long_name = "Nitrogen Mass in Model Ammonium Pool" ;
                nh4_mass:units = "g N/m2" ;
        float no3_mass(tile, iccp2, lat, lon) ;
                no3_mass:_FillValue = NaNf ;
                no3_mass:long_name = "Nitrogen Mass in Model Nitrate Pool" ;
                no3_mass:units = "g N/m2" ;
        float bnf_nat(tile, icc, lat, lon) ;
                bnf_nat:_FillValue = -999.f ;
                bnf_nat:long_name = "Natural Biological Nitrogen Fixation" ;
                bnf_nat:units = " kg N m$^{-2}$ s$^{-1}$" ;
        float bnf_ant(tile, icc, lat, lon) ;
                bnf_ant:_FillValue = -999.f ;
                bnf_ant:long_name = "Anthropogenic Biological Nitrogen Fixation" ;
                bnf_ant:units = " kg N m$^{-2}$ s$^{-1}$" ;
