%% Example script to implement functions in SolarRadiationModule
% Bryce Mihalevich
% Last modified: 5/3/18
%
% Description:

% Input Variables: 
%  xyPoints = [nx2] vector of query point(s) (x,y) in either UTM or
%  decimal degrees
%  xyCoord = string that equals 'UTM' or 'DD' to indicate coodinate
%  system of xyPoints
%  DEM_data = raster grid data (ideally a 10 km buffer around points of
%  interest; consider changing to greater or less buffer for certain areas)
%  DEM_ref = metadata for raster
%  azimuthDeg = [1xm] azimuth values to run the model for
%  h = number of meters to increase the ground surface elevation point of
%  interest. Removes some noise created by nearby cells but too much and
%  the illumination angles will be lower than they actually are. Use 0-3m
%  r = search radius in meters; 5000 seems to be sufficient for deep
%  canyons, 10000-15000 for open basin.
%  n = number of samples in search radius (if r/n = 1 == every 1 meter)
%  DOYs = days of the year to evaluate;
%  LocalTimes = times to evaluate (local time);
%  UTCoffset = timezone correction (Utah = -7)

clear; clc;
%% Change file paths here:
% Lat/lon points file (can be decimal degrees or UTM)
CORiv_latlon = '/Users/bryce/Documents/MATLAB/Colorado_River_Basin/GIS_data/SolarRadiationModule/Examples/GlenCan_CO_River_pnts.csv'; %lat long values in dec. deg.

% DEM data file with projected coordinates
DEMfilepath = '/Users/bryce/Documents/MATLAB/Colorado_River_Basin/GIS_data/SolarRadiationModule/Examples/GlenCan_dem_clip_UTM.tif';

%% Initialization Variables:
% ------------------ azimuth360.m function --------------------------
import csv
fid = fopen(CORiv_latlon,'r');
data = textscan(fid,'%f%f%f%[^\n\r]','Delimiter',',','HeaderLines',2); %read csv
fclose(fid);

riv_km = data{1}; %array of river km (not used in function)
lat_deg = data{2}; %array of river latitudes corresponding with river km
lon_deg = data{3}; %arrray of river longitudes corresponding with river km
xyCoord = 'DD'; %decimal degrees. If using utm coordinates change to "UTM"

% Example coordinates for UTM and DD - if not importing csv comment out lines above
% lat_deg = 36.93642; lon_deg = -111.48363; xyCoord = 'DD'; %Glen Canyon Dam
% utm_x = 456932.558709032; utm_y = 4087928.50129882; xyCoord = 'UTM'; %Glen Canyon Dam

h = 1; %meters
r = 10000;
n = 10000;

% import DEM
[DEM_data, DEM_ref] = geotiffread(DEMfilepath); %matrix of elevation data and reference/metadata 
azimuthDeg = 0:359; %degrees 


% ------------------ SolarInsolation.m function --------------------------
DOYs = [80,172,264,355]; %days of the year to evaluate;
LocalTimes =  0:(1/60):(23+56/60); %times during day to evaluate (local time);
UTCoffset = -7; %timezone correction

%% Run Illumination angle calculation function
[illum_angle,illum_elev,illum_dist]...
    = azimuth360([lon_deg,lat_deg],xyCoord,azimuthDeg,DEM_data,DEM_ref,h,r,n);

%% Run Solar Incidence function
[incidence, dateTimes, skyVF]...
    = SolarIncidence([lon_deg,lat_deg],xyCoord,azimuthDeg,illum_angle,DOYs,LocalTimes,UTCoffset);
