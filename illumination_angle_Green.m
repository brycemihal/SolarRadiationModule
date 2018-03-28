% Script for illumination angles for points within DEM
clear; clc; close all;


%% Import data
% Lat/lon points file
greenRiv_latlon = '/Users/bryce/Documents/School/USU Documents/PhD/Data/GIS Data/River_miles/Green_River_100m_pnts.csv'; %lat long values in dec. deg.
fid = fopen(greenRiv_latlon,'r');
data = textscan(fid,'%f%f%f%[^\n\r]','Delimiter',',','HeaderLines',2); %read csv
fclose(fid);

riv_km = data{1}; %array of river km
lat_deg = data{2}; %array of river latitudes corresponding with river km
lon_deg = data{3}; %arrray of river longitudes corresponding with river km

% Import DEM data file 
filepath = '/Users/bryce/Documents/GIS/ArcGIS/GreenRiver/New Folder/dem10m_green_10km.tif';
[DEM_data, DEM_ref] = geotiffread(filepath);

%% Modify DEM data if needed
% DEM data can be flipped. Use reference info to flip data
% if RowsStartFrom = 'east' then the data needs to be flipped left to right -- fliplr(DEM_data)
% if ColumnStartFrom = 'north' then the data is fliped upside down -- flipud(DEM_data)
if DEM_ref.ColumnsStartFrom == 'north'
    DEM_data = flipud(DEM_data);
end
if DEM_ref.RowsStartFrom == 'east'
    DEM_data = fliplr(DEM_data);
end

%% plot dem to check correct orientation
% DEM_data(DEM_data == -9999) = NaN;
% contourf(DEM_data(1:10000,1:5000))

%% Coordinate location(s) of query point(s)
% dem data coordinates are referenced in UTM, Point data is in lat lon
% Convert lat/lon dec. deg. to UTM 
[utm_x,utm_y,zone,isnorth] = utmups_fwd(lat_deg,lon_deg); %Array of utm x,y coordinates corresponding with river km
azimuthDeg = 0:359; % consider adding more degree increments

% utm_x = utm_x(1:1000); %use to run specific index values
% utm_y = utm_y(1:1000); %use to run specific index values

%% Run fuction to calculate illumination angles for given point or points
% illumination angle is a nx360 array where n = number of coordinate locations queried
[illum_angle,illum_elev,illum_dist,utm_elev,x_ind,y_ind,elev_ind] = azimuth360_v2([utm_x,utm_y],DEM_data,DEM_ref,azimuthDeg);

%% Test plot of dem (from single point)

% elev_ind(elev_ind == -9999) = NaN;
% xvals = squeeze(x_ind(end,:,:));
% yvals = squeeze(y_ind(end,:,:));
% elev_vals = squeeze(elev_ind(end,:,:));
% dem = contourf(xvals,yvals,elev_vals,10);
% 
% azimuthDeg = 0:359;
% azAngles = [0 90 180 270 360 450];
% radDeg = [90 0 -90 -180 -270 -360];
% azimuthDeg = round(interp1(azAngles,radDeg,azimuthDeg));
% 
% for i = 1:size(illum_dist,1)
%     xq = illum_dist(i,:).*cosd(azimuthDeg)+utm_x(i); %query point x
%     yq = illum_dist(i,:).*sind(azimuthDeg)+utm_y(i);
%     hold on
%     plot(xq,yq,'r')
%     plot(utm_x(i),utm_y(i),'ro')
% end


%% Test plot of illumination angles similar to yard plot
illum_angle2 = 90 - illum_angle';
figure(2)
contourf(riv_km,azimuthDeg,illum_angle2,100);
ax = gca;
ax.Children.EdgeColor = 'none';

%% custom color bar

% number of colors
% nc = ceil(max(max(illum_angle2)));
% 
% % 0 = white; mid = green; end = black
% c_map = [1,1,1 %white
%         1,1,0 %yellow
%         0,1,0 %green
%         0,0,1 %blue
%         0,0,0]; %black
% 
% n_yellow = round(nc/4);
% n_green = round(nc/2);
% n_blue = round(nc/2+nc/4);
% 
% 
% % white to green colorbar
% c_map = [linspace(1,0,100); linspace(1,0,100); linspace(1,0,100)]'


% colormap(flipud(pink)) 
%%
ax = gca;
ax.YLabel.String = 'Azimuth Angle';
ax.XLabel.String = 'River Kilometer';
c = colorbar;
c.Label.String = 'Illumination Angle';
c.Ticks = 10:10:80;
c.TickLabels = strcat(c.TickLabels,char(176));

ax.YTick = [0 90 180 270 359];
ax.YTickLabel(end) = {'360'};
ax.YTickLabel = strcat(ax.YTickLabel,char(176));

ax.TickDir = 'out';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.FontSize = 14;

ax.Title.String = 'Green River';
ax.Title.Position(1) = 675;

%%
t = text();
t.String = 'Gates of Lodore';
t.Position = [100 365 0];
t.HorizontalAlignment = 'center';
t.VerticalAlignment = 'bottom';
t.FontSize = 13;

t = text();
t.String = 'Desolation Canyon';
t.Position = [375 365 0];
t.HorizontalAlignment = 'center';
t.VerticalAlignment = 'bottom';
t.FontSize = 13;

t = text();
t.String = 'Labyrinth/Stillwater Canyon';
t.Position = [575 365 0];
t.HorizontalAlignment = 'center';
t.VerticalAlignment = 'bottom';
t.FontSize = 13;
