Instructions for calculating combined flood-heatwave events <br />
Heatwaves - based on dry bulb and wet bulb heatwaves. Calculated using Russo 2015 method. <br />
Floods - calcualated based on WAP method (Lu 2009 and Chen 2021) <br />

Inputs: CP4A and P25 model data, netcdf files (should work on any climate model data in netcdf format) <br />

Required variables: <br />
-temperature (2m) <br />
-pressure (surface) <br />
-specific humidity (surface) <br />
-preciptiation <br />

Process <br />
1. 	If not already, calcualte daily means from hourly data. <br />
	Ensure everything on same grid. <br />

	Scripts: <br />
	tas_prepare.py <br />
	precip_prepare.py <br />
	humidity_prepare.py <br />
	pressure_prepare.py <br />

2. 	Calculate wet bulb temperature.  <br />
	Currently using matlab version of Buzan's HumanIndexMod, with minor alterations to address dry areas not converging <br />

	Scripts (in wetbulb_matlab folder): <br />
	WetBulb.m (matlab version of HumanIndexMod, developed by Knopp) <br />
	run_wetBulb.m <br />

3. 	Calculate heatwaves <br />
	Based on Russo 2015 method <br />
	Inputs: <br />
	Td or Twb <br />

	Scripts (heatwaves): <br />
	calc_heatwaves.py <br />


4. 	Calculate WAP <br />
	Calculate weighted area precipitation, based on Lu 2009 method. <br />
	Currently using alpha = 0.9 and window size = 44 days, but script allows for change. <br />
	Alpha = 0.9, weight given to previous days. Amounts to precip after day 44 having very minimal weight (so don't bother including) <br />

	Inputs: <br />
	Daily pr <br />

	Scripts: <br />
	calc_wap3.py (wap folder) <br />

5.	Calculate floods <br />
	Calculate floods based on 95th percentile of WAP, based on Chen 2021 <br />

	Inputs:<br />
	WAP<br />

	Scripts (floods folder): <br />
	calc_floods.py<br />
	calc_floods_rcp85.py<br />

6.	Calculate combined events <br />
	Calculate combined flood-heatwave events, including sameday events, floods before heatwaves and floods after heatwaves <br />
	Currently calculating 7, 5, 3 days before/after heatwave, based on Chen 2021, but script allows this to change <br />

	Inputs (all inputs output from flood and heatwave scripts):<br />
	Flood days <br />
	Start day of flood events <br />
	End day of flood events <br />
	
	Heatwave days <br />
	Start of heatwave <br />
	End of heatwave <br />

	Scripts (in combined_events repo): <br />
	combiend events.py <br />

References: <br />

Chen, Y., Z. Liao, Y. Shi, Y. Tian, and P. Zhai. 2021. Detectable Increases in Sequential Flood-Heatwave Events Across China During 1961–2018. Geophysical Research Letters 48:1–10. <br />

Lu, E. 2009. Determining the start, duration, and strength of flood and drought with daily precipitation: Rationale. Geophysical Research Letters 36:1–5. <br />

Russo, S., J. Sillmann, and E. M. Fischer. 2015. Top ten European heatwaves since 1950 and their occurrence in the coming decades. Environmental Research Letters 10:124003. <br />


