%% Import data

%change to folder where WetBulb.m is (function to calculate wetbulb
%temperature)
cd /nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/floods_heatwaves/wetbulb_matlab

%% extract data to west africa

min_lat = 3.5;
max_lat = 20.0;
min_lon = -20.0;
max_lon = 16.0;


lons = ncread('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_historical.nc', 'longitude');
lats = ncread('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_historical.nc', 'latitude');
time = ncread('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_historical.nc', 'time');

timestart = 1;
timeend = 360;

latstart = find(lats > min_lat, 1, 'first');
latend = find(lats < max_lat, 1, 'last');

lonstart = find(lons > min_lon, 1, 'first');
lonend = find(lons < max_lon, 1, 'last');

start_loc = [lonstart latstart timestart];
end_loc = [(lonend-lonstart)+1 (latend-latstart)+1 timeend];

len_time = size(time);
len_time = len_time(1);

%% Data in order lons, lats, time

temp = ncread('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_p25_historical.nc', 'd03236', start_loc, end_loc);
pa = ncread('/nfs/a321/earsch/floods_heatwaves/input_data/pressure/pa_p25_daily_histo.nc', 'c00409', start_loc, end_loc);
hus = ncread('/nfs/a321/earsch/floods_heatwaves/input_data/humidity/spechumid_p25_daily_histo.nc', 'c03237', start_loc, end_loc);


%% Get into correct units for wetbulb func
% temp -> celsius
% pressure -> pa (already correct units)
% humidity -> kg/kg (already correct units)

temp = temp - 273.15;

%% Calculate wetbulb temp

disp('Computing wb')

[wb,Teq,epott] = WetBulb(temp, pa, hus,0, 1);

%% create dimensions for new netcdf file

disp('Creating file')

file_name = '/nfs/a321/earsch/floods_heatwaves/processed/wetbulb_temp/wb_p25_historical_year1.nc';
cmode = 'NETCDF4';

new_nc = netcdf.create(file_name, cmode);

%create dimensions
londim = netcdf.defDim(new_nc, 'longitude', end_loc(1));
latdim = netcdf.defDim(new_nc, 'latitude', end_loc(2));
timedim = netcdf.defDim(new_nc, 'time', timeend);

%create wb variable
wb_var = netcdf.defVar(new_nc, 'wb', 'NC_FLOAT', [londim, latdim, timedim]);
lon_var = netcdf.defVar(new_nc, 'longitude', 'NC_FLOAT', londim);
lat_var = netcdf.defVar(new_nc, 'latitude', 'NC_FLOAT', latdim);
time_var = netcdf.defVar(new_nc, 'time', 'NC_FLOAT', timedim);

%end define mode and enter data mode
netcdf.endDef(new_nc);

%add data to netcdf file
new_lons = lons(lonstart:lonend);
new_lats = lats(latstart:latend);
new_times = time(timestart:timeend);

netcdf.putVar(new_nc, wb_var, wb);
netcdf.putVar(new_nc, lon_var, new_lons);
netcdf.putVar(new_nc, lat_var, new_lats);
netcdf.putVar(new_nc, time_var, new_times);

netcdf.close(new_nc);

