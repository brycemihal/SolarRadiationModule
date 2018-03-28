%% Function to calculate 0:360 degree asimuth angles for a lat/lon point
function [illum_angle,illum_elev,illum_dist,utm_elev,x_ind,y_ind,elev_ind] = azimuth360(xyPoints,DEM_data,DEM_ref,azimuthDeg)

% Input Variables
%       latlonPoint = [nx2] vector of query point(s) (x,y) in UTM
%       DEM_data = raster grid data (10 km buffer around points of interest; consider changing to greater or less buffer for certain areas)
%       DEM_ref = metadata for raster
%       azimuthDeg = [1xm] vector of values

% Output Variables
%       illum_angle = [n x m] vector of illumination agles where n = number
%       of latlongPoints (from vertical)
%       utm_elev = [n x 1] vector with elevations of query points (meter)
%       illum_elev = [n x m] elevation at the point with greatest illumination angle (meters)
%       illum_dist = [n x m] distance from query point to the greatest illumination angle (meters)
%       x_ind = [n x m x qdist] matrix of all queired x points
%       x_ind = [n x m x qdist] matrix of all queired y points
%       elev_ind = [n x m x qdist] matrix of all queired elevations

%% Calculate the illumination angle for each lat/lon point

% set search radius; consider changing this to query every intersected cell (https://www.mathworks.com/matlabcentral/answers/230155-how-to-determine-which-grid-cells-a-line-segment-passes-through)
r = 10000; %search radius in meters
n = 10000; % number of samples in search radius (r/n = 1 == every 1 meter)
qdist = linspace(0,r,n+1);%query distances

% matlab has 0 degrees == 0*pi (east) and increases counterclockwise
% azimuth angles start at 0 degrees = north and increase counterclockwise
% convert angles using interp
azAngles = [0 90 180 270 360 450];
radDeg = [90 0 -90 -180 -270 -360];
azimuthDegq = round(interp1(azAngles,radDeg,azimuthDeg),1);

% Raster cell size
xCellSize = DEM_ref.CellExtentInWorldX;
yCellSize = DEM_ref.CellExtentInWorldY;

for pnt = 1:size(xyPoints,1) % loop for each point in array
    %% Elevation of query point
    % Determine elevation from DEM grid        
    col = ceil(abs((DEM_ref.XWorldLimits(1)-xyPoints(pnt,1))/xCellSize));
    row = ceil(abs((DEM_ref.YWorldLimits(1)-xyPoints(pnt,2))/yCellSize));
    elev_m = DEM_data(row,col);
    
    % interpolate row/column index (can this be done without interpolation/faster?)
%     col = ceil(interp1(DEM_ref.XWorldLimits,DEM_ref.XIntrinsicLimits,xyPoints(pnt,1))-.5); %(column)
%     row = ceil(interp1(DEM_ref.YWorldLimits,DEM_ref.YIntrinsicLimits,xyPoints(pnt,2))-.5); %(row)
%     elev_m = DEM_data(row,col); %elevation in meters of x/y
    
    % loop for azimuth angles and determine cell value for each distance
    for i=1:length(azimuthDegq) % number of azimuth angles
        
        xq = qdist.*cosd(azimuthDegq(i))+xyPoints(pnt,1); %query point x
        yq = qdist.*sind(azimuthDegq(i))+xyPoints(pnt,2); %query point y
        
        % interpolate row/column index
%         colq = ceil(interp1(DEM_ref.XWorldLimits,DEM_ref.XIntrinsicLimits,xq)-.5); %(column)
%         rowq = ceil(interp1(DEM_ref.YWorldLimits,DEM_ref.YIntrinsicLimits,yq)-.5); %(row)

        colq = ceil(abs((DEM_ref.XWorldLimits(1)-xq)/xCellSize)); % gets the column index for x coordination location
        rowq = ceil(abs((DEM_ref.YWorldLimits(1)-yq)/yCellSize)); % gets the row index for y coordination location
        
        % get elevation values for each query point
%         for j = 1:length(qdist)
%                 elev_mq(pnt,j) = DEM_data(row(j),col(j));%queried elvation in meters at new xy coordinates
%         end
        % use linear indexing to assign DEM_data elevation to array (faster than for loop?)
        elev_mq(pnt,:) = DEM_data(sub2ind(size(DEM_data),rowq(:),colq(:))); %queried elevation in meters at new xy coordinates
        
        % calculate illumination anlge
        dH = elev_mq(pnt,:)-elev_m; %change in height
        [max_illum_a,loc] = max(atand(dH(2:end)./qdist(2:end))); % max angle
        illum_angle(pnt,i) = 90-max_illum_a; % subtract 90 so angle is from the vertical
        
        % additional metrics
        utm_elev(pnt) = elev_m;
        x_ind(pnt,i,:) = xq; %x coord for each dist %(Slows down code, comment out if not needed)
        y_ind(pnt,i,:) = yq; %y coord for each dist %(Slows down code, comment out if not needed)
        elev_ind(pnt,i,:) = elev_mq(pnt,:); % elevation for each distance %(Slows down code, comment out if not needed)
        illum_elev(pnt,i) = elev_mq(loc); % rim elevation for max illumination angle
        illum_dist(pnt,i) = qdist(loc); % distance to rim elevation for max illum angle
    end
end

% utm_elev = [];
% x_ind = [];
% y_ind = [];
% elev_ind = [];



