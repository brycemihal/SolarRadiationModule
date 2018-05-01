%% Compare Lat lon points from Yard
clear; clc; close all;
% Lat/lon points file
CORiv_latlon = '/Users/bryce/Documents/School/USU Documents/PhD/Data/GIS Data/River_miles/CO_River_Yard_pnts.csv'; %lat long values in dec. deg.
fid = fopen(CORiv_latlon,'r');
data = textscan(fid,'%f%f%f%[^\n\r]','Delimiter',',','HeaderLines',2); %read csv
fclose(fid);

riv_km = data{1}; %array of river km
lat_deg = data{2}; %array of river latitudes corresponding with river km
lon_deg = data{3}; %arrray of river longitudes corresponding with river km
xyCoord = 'DD'; %decimal degrees
h = 2; %meters

% Import DEM data file with projected coordinates
filepath = '/Users/bryce/Documents/MATLAB/Colorado_River_Basin/GIS_data/SolarRadiationModule/Examples/GlenCan_dem_clip_UTM.tif';
[DEM_data, DEM_ref] = geotiffread(filepath);

% Run function
azimuthDeg = 0:359;                                                             % (x)    (y)
[illum_angle,illum_elev,illum_dist,utm_elev,x_ind,y_ind,elev_ind] = azimuth360([lon_deg,lat_deg],xyCoord,h,DEM_data,DEM_ref,azimuthDeg,10000,10000,0);


%%
figure
DEM_data2 = flipud(DEM_data);
DEM_data2(DEM_data2<0) = NaN;
surf(DEM_data2)
ax = gca;
ax.Children.EdgeColor = 'none';
ax.View = [0 90];



%% import data from yard model

yard_file = '/Users/bryce/Documents/School/USU Documents/PhD/Data/GIS Data/LightModel/yard_illum_angle.csv';
fid = fopen(yard_file,'r');
formatSpec = strcat(repmat('%f',1,201),'%[^\n\r]');
yard_data = textscan(fid,formatSpec,'Delimiter',',','HeaderLines',2);

% convert cells to array
yard_illum_angle = cell2mat(yard_data(2:end-1))'; % row == riv km
yard_azimuthDeg = yard_data{1}';

%% calculate statistics
NSE = @(x1,x2) 1-(sum((x2-x1).^2)/sum((x1-mean(x1)).^2));  % Nash - Sutcliffe medel efficiency coefficient
RMSE = @(x1,x2) sqrt(mean((x1-x2).^2)); %root mean square error
for i = 1:10:size(illum_angle,1)
    rmse_calc(i) = RMSE(yard_illum_angle(i,:),illum_angle(i,:));
end
%% plot comparison of models

for i = 1:10:size(illum_angle,1)
    figure(i)
    hold on
    yard_h = plot(yard_azimuthDeg,yard_illum_angle(i,:),'r');
    my_h = plot(azimuthDeg,illum_angle(i,:),'b');
    hh = gobjects(1);
    ax = gca;
    ax.XLim = [0 360];
    ax.XLabel.String = 'Azimuth Degrees';
    ax.YLabel.String = 'Illumination angle (from Horizontal)';
    ax.Title.String = sprintf('CO River %2.1f km below Glen Canyon Dam',riv_km(i));
    ax.FontSize = 13;
    l = legend;
    l.String = {'Yard model','My model','RMSE = .2'};
    l.Location = 'southeast';
    l.FontSize = 11;
    t = text;
    t.String = sprintf('RMSE = %2.3f',rmse_calc(i));
    t.Position = [max(ax.XLim)-3 max(ax.YLim) 0];
    t.HorizontalAlignment = 'right';
    t.VerticalAlignment = 'top';
    t.FontSize = 11;
    
    save_eps = sprintf('/Users/bryce/Documents/MATLAB/Colorado_River_Basin/GIS_data/Illum_angle_comparison/illum_angle_comp_km%2.1f.eps',riv_km(i));
    save_tif = sprintf('/Users/bryce/Documents/MATLAB/Colorado_River_Basin/GIS_data/Illum_angle_comparison/illum_angle_comp_km%2.1f.tif',riv_km(i));
    saveas(gcf,save_eps,'epsc')
    saveas(gcf,save_tif,'tiffn')
    close all
end
