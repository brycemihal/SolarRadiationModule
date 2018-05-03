%% Function to calculate 0:360 degree asimuth angles for a lat/lon points

function [illum_angle,illum_elev,illum_dist,utm_elev,x_ind,y_ind,elev_ind] = azimuth360(xyPoints,xyCoord,h,DEM_data,DEM_ref,azimuthDeg,r,n,return_grid)
% Bryce Mihalevich
% Last modified: 3/4/18
%
%
% Description: 
% This function will calculate the illumination angle (the angle from the
% horizontal between a point on a surface to the highest topographic point 
% within line of site) given an input xy coordinate point, DEM raster and 
% azimuth angles. xy coordinate points can be either UTM or decimal
% degrees but must be specified in the xyCoord variable. The coordinate 
% system of the DEM must be a projected coordinate system such as NAD83 
% zone 12.
% 
% Input Variables
%       xyPoints = [nx2] vector of query point(s) (x,y) in either UTM or
%       decimal degrees
%       xyCoord = string that equals 'UTM' or 'DD' to indicate coodinate
%       system
%       h = meters above ground surface. Model seems to work better with
%       observation point 1-3 meters above surface.
%       DEM_data = raster grid data (10 km buffer around points of interest
%       consider changing to greater or less buffer for certain areas)
%       DEM_ref = metadata for raster
%       azimuthDeg = [1xm] vector of values
%       r = search radius in meters; 5000 seems to be sufficient for deep
%       canyons, 10000-15000 for open basin.
%       n = number of samples in search radius (r/n = 1 == every 1 meter)
%       return_grid = binary 1 or 0. If 1 additional data will be available
%       in the output. Good making plots of raster. Slows down code a lot.
%       Only recommend if doing < 10 points.
% 
% Output Variables
%       illum_angle = [n x m] vector of illumination agles where n = number
%       of latlongPoints (measured from the HORIZONTAL)
%       utm_elev = [n x 1] vector with elevations of query points (meter)
%       illum_elev = [n x m] elevation at the point with greatest 
%       illumination angle (meters)
%       illum_dist = [n x m] distance from query point to the greatest 
%       illumination angle (meters)
%       x_ind = [n x m x qdist] matrix of all queired x points
%       x_ind = [n x m x qdist] matrix of all queired y points
%       elev_ind = [n x m x qdist] matrix of all queired elevations
% 
% %% Example 1
% %Lat/lon points
% lat_deg = 36.860720;
% lon_deg = -111.608441;
% xyCoord = 'DD'; %decimal degrees
% h = 2; %meters
% 
% %Import DEM data file with projected coordinates
% filepath = '/Users/bryce/Documents/MATLAB/Colorado_River_Basin/GIS_data/SolarRadiationModule/Examples/GlenCan_dem_clip_UTM.tif';
% [DEM_data, DEM_ref] = geotiffread(filepath);
% 
% azimuthDeg = 0:359;                                                             % (x)    (y)
% [illum_angle,illum_elev,illum_dist,utm_elev,x_ind,y_ind,elev_ind] = azimuth360([lon_deg,lat_deg],xyCoord,h,DEM_data,DEM_ref,azimuthDeg,10000,10000,0);

%% Flip DEM data if needed
% code is writen so y data starts in the south; DEM's typically start in the north. 
% DEM data can be flipped. Use reference info to flip data
% if RowsStartFrom = 'east' then the data needs to be flipped left to right -- fliplr(DEM_data)
% if ColumnStartFrom = 'north' then the data is fliped upside down -- flipud(DEM_data)
if strcmp(DEM_ref.ColumnsStartFrom,'north')
    DEM_data = flipud(DEM_data);
    DEM_ref.ColumnsStartFrom = 'south';
end
if strcmp(DEM_ref.RowsStartFrom,'east')
    DEM_data = fliplr(DEM_data);
    DEM_ref.RowsStartFrom = 'west';
end

%% Determine coordinate system of points and convert xyCoords if needed
% if xy coordinates do not match raster, convert xy to match
if strcmp(xyCoord,'DD') %then convert to UTM
        [utm_x,utm_y] = utmups_fwd(xyPoints(:,2),xyPoints(:,1)); %DD to UTM conversion
        xyPoints = [utm_x,utm_y];
end

% Raster cell size
xCellSize = DEM_ref.CellExtentInWorldX;
yCellSize = DEM_ref.CellExtentInWorldY;

% limits of DEM
xLimits = DEM_ref.XWorldLimits;
yLimits = DEM_ref.YWorldLimits;
    
%% Calculate the illumination angle for each lat/lon point
% set search radius; consider changing this to query every intersected cell (https://www.mathworks.com/matlabcentral/answers/230155-how-to-determine-which-grid-cells-a-line-segment-passes-through)
qdist = linspace(0,r,n+1); %query distances in meters

% matlab has 0 degrees == 0*pi (east) and increases counterclockwise
% azimuth angles start at 0 degrees = north and increase counterclockwise
% convert angles using interp
azAngles = [0 90 180 270 360 450];
radDeg = [90 0 -90 -180 -270 -360];
azimuthDegq = round(interp1(azAngles,radDeg,azimuthDeg),1);

%% Preallocate memory to variables for speed
% elev_mq
% illum_angle
% illum_elev
% illum_dist

%% Loop for number of xy points
for pnt = 1:size(xyPoints,1) % loop for each point in array
    %% Elevation of query point
    % Determine elevation from DEM grid        
    col = ceil(abs((xLimits(1)-xyPoints(pnt,1))/xCellSize));
    row = ceil(abs((yLimits(1)-xyPoints(pnt,2))/yCellSize));
    elev_m = DEM_data(row,col)+h; %the model seems to do better with the point 1-3 meters above he DEM surface
    
    % loop for azimuth angles and determine cell value for each distance
    for i=1:length(azimuthDegq) % number of azimuth angles
        
        xq = qdist.*cosd(azimuthDegq(i))+xyPoints(pnt,1); %query point x (easting in meters)
        yq = qdist.*sind(azimuthDegq(i))+xyPoints(pnt,2); %query point y (northing in meters)

        colq = ceil(abs((xLimits(1)-xq)/xCellSize)); % gets the column index for x (easting) locations
        rowq = ceil(abs((yLimits(1)-yq)/yCellSize)); % gets the row index for y (northing) locations
        
        % use linear indexing to assign DEM_data elevation to array (faster than for loop?)
        elev_mq(pnt,:) = DEM_data(sub2ind(size(DEM_data),rowq(:),colq(:))); %queried elevation in meters at new xy coordinates (this line breaks the code sometimes...)

        % calculate illumination anlge
        dH = elev_mq(pnt,:)-elev_m; %change in height
        [illum_angle(pnt,i),loc] = max(atand(dH(2:end)./qdist(2:end))); % max angle
        illum_elev(pnt,i) = elev_mq(loc); % rim elevation for max illumination angle
        illum_dist(pnt,i) = qdist(loc); % distance to rim elevation for max illum angle
        
        % additional metrics (slows down code significantly)
        if return_grid == 1
            utm_elev(pnt) = elev_m;
            x_ind(pnt,i,:) = xq; %x coord for each dist %(Slows down code, comment out if not needed)
            y_ind(pnt,i,:) = yq; %y coord for each dist %(Slows down code, comment out if not needed)
            elev_ind(pnt,i,:) = elev_mq(pnt,:); % elevation for each distance %(Slows down code, comment out if not needed)
        end
    end
end

% if additional metrics are not needed set grid_data to 0 in function
if return_grid == 0
    utm_elev = [];
    x_ind = [];
    y_ind = [];
    elev_ind = [];
end
