%% Compare model to Yard results
clear; clc; close all;

%% Initialize variables
% % Lat/lon points file
CORiv_latlon = '/Users/bryce/Documents/School/USU Documents/PhD/Data/GIS Data/River_miles/CO_River_Yard_pnts.csv'; %lat long values in dec. deg.
fid = fopen(CORiv_latlon,'r');
data = textscan(fid,'%f%f%f%[^\n\r]','Delimiter',',','HeaderLines',2); %read csv
fclose(fid);

riv_km = data{1}; %array of river km
lat_deg = data{2}; %array of river latitudes corresponding with river km
lon_deg = data{3}; %arrray of river longitudes corresponding with river km
% lat_deg = 36.92807; lon_deg = -111.47948; %1km below Glen Canyon Dam
% lat_deg = 36.93642; lon_deg = -111.48363; %Glen Canyon Dam
xyCoord = 'DD'; %decimal degrees
h = 1; %meters

% Import DEM data file with projected coordinates
filepath = '/Users/bryce/Documents/MATLAB/Colorado_River_Basin/GIS_data/SolarRadiationModule/Examples/GlenCan_dem_clip_UTM.tif';
[DEM_data, DEM_ref] = geotiffread(filepath);

%% Run Illumination angle calculation function
azimuthDeg = 0:359;                                        % (x)    (y)
[illum_angle,illum_elev,illum_dist,~,~,~,~] = azimuth360([lon_deg,lat_deg],xyCoord,h,DEM_data,DEM_ref,azimuthDeg,10000,10000,0);

%% Initialize variables

% Day of year to test (julian days): spring equinox (80), 
% summer solstice (172), fall equinox (264), winter solstice (355) 
DOYs = [80,172,264,355]; %days of the year to evaluate;
LocalTimes =  0:(1/60):(23+56/60); %times during day to evaluate (local time);
UTCoffset = -7; %timezone correction

%% Run Solar Incidence function
[sunrise,sunriseAz,sunset,sunsetAz,TotalDayLength,incidence,directBeamDayLength] = SolarIncidence([lon_deg,lat_deg],xyCoord,azimuthDeg,illum_angle,DOYs,LocalTimes,UTCoffset);

