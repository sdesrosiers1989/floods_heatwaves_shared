# Instructions for calculating combined flood-heatwave events <br />
Heatwaves - based on dry bulb and wet bulb heatwaves. Calculated using Russo 2015 method. <br />
Hotdays - based on 95th percentile of daily average temperature. Alternative to using heatwaves. <br />
Floods - calcualated based on WAP method (Lu 2009 and Chen 2021). <br />

Inputs:
- CP4A and P25 model data, netcdf files (should work on any climate model data in netcdf format) 
- Required variables: <br />
  - temperature (2m) 
  - pressure (surface) 
  - specific or relative humidity (surface) 
  - preciptiation 

## Process <br />
1. Calculate daily means from hourly data. <br />
	Ensure everything on same grid. <br />

	Scripts (note: these are within the private floods_heatwave repos - email for access): <br />
	tas_prepare.py <br />
	precip_prepare.py <br />
	humidity_prepare.py <br />
	pressure_prepare.py <br />

2. Calculate wet bulb temperature.  <br />
	Currently using matlab version of Buzan's HumanIndexMod, with minor alterations to address dry areas not converging <br />

	Scripts (in wetbulb_matlab folder): <br />
	WetBulb.m (matlab version of HumanIndexMod, developed by Knopp) <br />
	run_wetBulb.m - run this to apply the WetBulb.m function to climate data <br />
	wetbulb_calc_concatmatlabfiles.py - if run_wetBulb.m was applied to smaller sections of climate data, run this to put the files back together into one file <br />

3. Calculate heatwaves/hotdays <br />
	Heatwaves based on Russo 2015 method <br />
	
	Inputs: <br />
	Td or Twb <br />

	Scripts (heatwaves): <br />
	calc_heatwaves.py <br />
	
	Alternative: calculate hot days instead (note: unlike heatwaves, the historical threshold for hotdays is also used in the future) <br />
	Scripts (hotdays folder): <br />
	calc_hotdays.py <br />
	hotdays_panafrica_concatfiles.py - if hotdays was applied to smaller sections of climate data, run this to put the files back together into one file <br />


4. Calculate WAP <br />
	Calculate weighted area precipitation, based on Lu 2009 method. <br />
	Currently using alpha = 0.9 and window size = 44 days, but script allows for change. <br />
	Alpha = 0.9, weight given to previous days. Amounts to precip after day 44 having very minimal weight (so don't bother including) <br />
	Pan-africa: WAP calculated in sections and saved (otherwise takes too long) - files then concatenated into one <br />

	Inputs: <br />
	Daily pr <br />

	Scripts (wap folder): <br />
	calc_wap3.py <br />
	wap_panafrica_concatfiles.py <br />

5. Calculate floods <br />
	Calculate floods based on 95th percentile of WAP, and having WAP >= 10.0 (precentile threshold based on Chen 2021) <br />

	Inputs:<br />
	WAP<br />

	Scripts (floods folder): <br />
	flood_functions <br /> - calc_floods call this to calculate floods
	calc_floods.py <br />
	calc_floods_rcp85.py<br />
	floods_panafrica_concatfiles.py - put floods files back together <br/>

6. Calculate combined events <br />
	Calculate combined flood-heatwave/hotday events, including sameday events, floods before heatwaves and floods after heatwaves <br />
	Currently calculating 1 - 7, and 15 days before/after heatwave, but script allows this to change <br />

	Inputs (all inputs output from flood and heatwave/hotday scripts):<br />
	Flood days <br />
	Start day of flood events <br />
	End day of flood events <br />
	
	Heatwave days <br />
	Start of heatwave <br />
	End of heatwave <br />

	Scripts (in combined_events folder): <br />
	combined_events.py <br />
	
	Alternative: combined flood hotday events instead <br />
	Scripts (in combined_events folder): <br />
	combined_events_hotdays.py <br />
	combinedevents_panafrica_concatfiles.py - run this to put the files back together into one file <br />

	

## Progress
### West Africa (CP4A and P25)
- [x] Wetbulb temperature
- [x] Heatwaves
- [x] WAP
- [x] Floods
- [x] Combined events

### Pan Africa (CP4A and P25)
- [X] Wetbulb temperature
- [X] Heatwaves
- [x] Hot days
- [X] WAP
- [X] Floods
- [ ] Combined events

## Still to do
### Global (ERA5 and GPM)
- [ ] Relative humidity
- [ ] Wetbulb temperature
- [ ] Heatwaves
- [ ] Hot days
- [ ] WAP
- [ ] Floods
- [ ] Combined events

## References <br />

Chen, Y., Z. Liao, Y. Shi, Y. Tian, and P. Zhai. 2021. Detectable Increases in Sequential Flood-Heatwave Events Across China During 1961–2018. Geophysical Research Letters 48:1–10. <br />

Lu, E. 2009. Determining the start, duration, and strength of flood and drought with daily precipitation: Rationale. Geophysical Research Letters 36:1–5. <br />

Russo, S., J. Sillmann, and E. M. Fischer. 2015. Top ten European heatwaves since 1950 and their occurrence in the coming decades. Environmental Research Letters 10:124003. <br />


