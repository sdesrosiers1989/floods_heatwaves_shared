Instructions for calculating combined flood-heatwave events < /b>
Heatwaves - based on dry bulb and wet bulb heatwaves. Calculated using Russo 2015 method.
Floods - calcualated based on WAP method (Lu 2009 and Chen 2021)

Inputs: CP4A and P25 model data, netcdf files (should work on any climate model data in netcdf format)

Required variables:
-temperature (2m)
-pressure (surface)
-specific humidity (surface)
-preciptiation

Process
1. 	If not already, calcualte daily means from hourly data.
	Ensure everything on same grid.

	Scripts:
	tas_prepare.py
	precip_prepare.py
	humidity_prepare.py
	pressure_prepare.py

2. 	Calculate wet bulb temperature. 
	Currently using matlab version of Buzan's HumanIndexMod, with minor alterations to address dry areas not converging

	Scripts (in wetbulb_matlab repo):
	WetBulb.m (matlab version of HumanIndexMod, developed by Knopp)
	run_wetBulb.m

3. 	Calculate heatwaves
	Based on Russo 2015 method
	Inputs:
	Td or Twb

	Scripts (heatwaves repo):
	calc_heatwaves.py


4. 	Calculate WAP
	Calculate weighted area precipitation, based on Lu 2009 method.
	Currently using alpha = 0.9 and window size = 44 days, but script allows for change.
	Alpha = 0.9, weight given to previous days. Amounts to precip after day 44 having very minimal weight (so don't bother including)

	Inputs:
	Daily pr

	Scripts:
	calc_wap3.py (wap repo)

5.	Calculate floods
	Calculate floods based on 95th percentile of WAP, based on Chen 2021

	Inputs:
	WAP

	Scripts (floods repo):
	calc_floods.py
	calc_floods_rcp85.py

6.	Calculate combined events
	Calculate combined flood-heatwave events, including sameday events, floods before heatwaves and floods after heatwaves
	Currently calculating 7, 5, 3 days before/after heatwave, based on Chen 2021, but script allows this to change

	Inputs (all inputs output from flood and heatwave scripts):
	Flood days
	Start day of flood events
	End day of flood events
	
	Heatwave days
	Start of heatwave
	End of heatwave

	Scripts (in combined_events repo):
	combiend events.py

References:

Chen, Y., Z. Liao, Y. Shi, Y. Tian, and P. Zhai. 2021. Detectable Increases in Sequential Flood-Heatwave Events Across China During 1961–2018. Geophysical Research Letters 48:1–10.

Lu, E. 2009. Determining the start, duration, and strength of flood and drought with daily precipitation: Rationale. Geophysical Research Letters 36:1–5.

Russo, S., J. Sillmann, and E. M. Fischer. 2015. Top ten European heatwaves since 1950 and their occurrence in the coming decades. Environmental Research Letters 10:124003.


