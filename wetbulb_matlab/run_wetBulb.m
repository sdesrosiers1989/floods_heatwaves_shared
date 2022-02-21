%% Calculate WetBulb temperature
% Apply WetBulb.m function to input model data (temp, pressure, specific
% humidity) and save output as a netcdf file

%% Set up functions
%change to folder where WetBulb.m is (function to calculate wetbulb
%temperature)
cd /nfs/see-fs-02_users/earsch/Documents/Leeds/Repos/floods_heatwaves/floods_heatwaves_shared/wetbulb_matlab;

%% Set up input data
% set model and scen, and then input file names

%area
%area = 'wa';
area = 'pa'; %pan-africa

%model
model = 'p25';
model = 'cp4';

%set year (note matlab starts at 1, indexing is inclusive)
%note that temp goes up to 3660 (feb 2007) while pa and RHs only go up to
%3600 (end up 2006) -> doesn't matter as drop last half of 2006 anyway, but
%be aware that time length is different
time_name = 'part1';
timestart = 1;
timeend = 3600;

%scenarios
scen = 'histo';
t_scen = 'historical';
scen = 'rcp85';
t_scen = 'rcp85';

%grid
if strcmp(model,'cp4') == true
    grid = '_p25grid';
else
    grid = '';
end

%files
temp_file = strcat('/nfs/a321/earsch/Tanga/Data/CP4_Processed/tas/tas_day_', model, '_', t_scen, grid, '.nc');
pa_file = strcat('/nfs/a321/earsch/floods_heatwaves/input_data/pressure/pa_', model, grid, '_daily_', scen, '.nc');
hus_file = strcat('/nfs/a321/earsch/floods_heatwaves/input_data/humidity/spechumid_', model, grid, '_daily_', scen, '.nc');

% save info
save_path = '/nfs/a321/earsch/floods_heatwaves/processed/wetbulb_temp/wb_';
file_name = strcat(save_path, model, '_', scen, '_', time_name, '_', area, '.nc');



%% extract data to west africa

min_lat = 3.5;
max_lat = 20.0;
min_lon = -20.0;
max_lon = 16.0;


lons = ncread(temp_file, 'longitude');
lats = ncread(temp_file, 'latitude');
time = ncread(temp_file, 'time');

len_lons =  size(lons)
len_lats = size(lats)

latstart = find(lats > min_lat, 1, 'first');
latend = find(lats < max_lat, 1, 'last');

lonstart = find(lons > min_lon, 1, 'first');
lonend = find(lons < max_lon, 1, 'last');

start_loc = [lonstart latstart timestart];
end_loc = [(lonend-lonstart)+1 (latend-latstart)+1 (timeend - timestart)+1];

len_time = size(time);
len_time = len_time(1);

%grid
if strcmp(area,'pa') == true
    latstart = 1
    lonstart = 1
    lonend = len_lons(1)
    latend = len_lats(1)
    
    start_loc = [1 1 timestart];
    end_loc = [Inf Inf (timeend - timestart) + 1];
end


%% Import data, extracting to west africa

temp = ncread(temp_file, 'd03236', start_loc, end_loc);
pa = ncread(pa_file, 'c00409', start_loc, end_loc);
hus = ncread(hus_file, 'c03237', start_loc, end_loc);


%% Get into correct units for wetbulb func
% temp -> celsius
% pressure -> pa (already correct units)
% humidity -> kg/kg (already correct units)

temp = temp - 273.15;

%% Calculate wetbulb temp

disp('Computing wb')

%for y = 1:103
%    for z = 1:70
%        disp_var = [y z];
%        disp(disp_var)
%        t = temp(y, z, 1);
%        p = pa(y, z, 1);
%        h = hus(y, z, 1);
%        [wb, Teq, epott] = WetBulb(t, p, h, 0, 1);
%    end
%end

[wb,Teq,epott] = WetBulb(temp, pa, hus, 0, 1);

%% create dimensions for new netcdf file

disp('Creating file')

cmode = 'NETCDF4'; %will overwrite existing files, create netcdf4 file

new_nc = netcdf.create(file_name, cmode);

%create dimensions
londim = netcdf.defDim(new_nc, 'longitude', lonend);
latdim = netcdf.defDim(new_nc, 'latitude', latend);
timedim = netcdf.defDim(new_nc, 'time', end_loc(3));

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

%%

netcdf.putVar(new_nc, wb_var, wb);
netcdf.putVar(new_nc, lon_var, new_lons);
netcdf.putVar(new_nc, lat_var, new_lats);
netcdf.putVar(new_nc, time_var, new_times);

netcdf.close(new_nc);

disp('Complete')

