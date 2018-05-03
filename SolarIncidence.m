function[sunriseTime,sunriseAz,sunsetTime,sunsetAz,TotalDayLength,incidence,directBeamDayLength] = SolarIncidence(xyPoints,xyCoord,azimuthDeg,illum_angle,DOYs,LocalTimes,UTCoffset)
% Bryce Mihalevich
% Last Modified: 5/2/18
% 
% NOTE: Function is incomplete; Model is working but output variables are
% not yet added to function output
% ADD example code
% FIX reshape incidence code
%
% Description: function to determine solar incidence
%
% Input Variables
%       xyPoints = [nx2] vector of query point(s) (x,y) in either UTM or
%       decimal degrees
%       xyCoord = string that equals 'UTM' or 'DD'
%       azimuthDeg = [1xm] vector of values
%       illum_angle = [n x m] vector of illumination agles where n = number
%       of latlongPoints (measured from the HORIZONTAL)
%       r = search radius in meters; 5000 seems to be sufficient for deep
%       canyons, 10000-15000 for open basin.
%       n = number of samples in search radius (r/n = 1 == every 1 meter)
%       return_grid = binary 1 or 0. If 1 additional data will be available
%       in the output. Good making plots of raster. Slows down code a lot.
%       Only recommend if doing < 10 points.
%       DOYs = [1xn] Day of Year (integers): Jan. 1 = 1, Dec. 31 = 365
%       Localtime: local standard time in decimal hours: 0.0=midnight, 23.0=11pm
%       UTCoffset (hours): Time shift relative to UTC (hours) (i.e. Utah, -7)
% 
% Output Variables

%       incidence
%       Celestial variables:
%       sunriseTime = calculated time of sunrise specific to a day and point
%       sunriseAz = azimuth angle at sunrise specific to a day and point
%       sunriseZe = zenith angle at sunrise specific to a day and point
%       sunsetTime = calculated time of sunset specific to a day and point
%       sunsetAz = azimuth angle at sunset specific to a day and point
%       sunsetZe = zenith angle at sunset specific to a day and point
%       TotalDayLength = sunrise minus sunset in minutes
%       
%       Topo variables:
%       TopoRiseTime = time when sun first hist point of interest  
%       TopoRiseSolAz = azimuth angle when sun first hits point of interst
%       TopoRiseSolAlt = 90-zenith angle when sun first hits point of interest 
%       TopoSetTime = time when sun first hist point of interest  
%       TopoSetSolAz = azimuth angle when sun first hits point of interst
%       TopoSetSolAlt = 90-zenith angle when sun first hits point of interest
%       directBeamDayLength = TopoSetTime minut TopoRiseTime in minutes

%       HemisphericalArea
%       AzimuthOfLowestAltitude90to269
%       AzimuthOfLowestAltitude270to89
%       CanyonAzimuth
%       CanyonHemiAngle

%       Radiation variabels:
%       diffuseRiseTime = time of morning civil twilight, defined as 6 degrees
%       (from the vertical) below the horizon (i.e.: ~96 degrees)
%       diffuseSetTime = time of evening civil twilight, defined as 6 degrees
%       (from the vertical) below the horizon  (i.e.: ~96 degrees)
%       difuseLength = diffusesetTime minus diffuseRiseTime in minutes
%       InsolationDirect = ?
%       InsolationDiffuse = ?
%       InsolationTotal = InsolationDiffuse plus InsolationDirect

% %% Example 1

%% Determine coordinate system of points and convert xyCoords if needed
% if xy coordinates are in UTM then convert
if strcmp(xyCoord,'UTM') %then convert to UTM
        zone = 12;
        isnorth = 1;
        [lat,lon] = utmups_inv(xyPoints(:,1),xyPoints(:,2),zone,isnorth); %UTM to DD conversion
        xyPoints = [lon,lat]; %CHECK THIS
end

% separate x and y values
lon_deg = xyPoints(:,1);
lat_deg = xyPoints(:,2);

%% preallocate memory to variables for speed
% Arrays
SkyVF = zeros(size(xyPoints,1),1); %floating point array
sunriseTime = cell(size(xyPoints,1),size(DOYs,2)); %character array
sunriseAz = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
sunriseZe = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
TotalDayLength = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
sunsetTime = cell(size(xyPoints,1),size(DOYs,2)); %character array
sunsetAz = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
sunsetZe = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
incidence = false(size(DOYs,2),size(LocalTimes,2),size(xyPoints,1)); %3D logical array [DOY x incidence x point]
directBeamDayLength = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoRiseTime = cell(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoRiseSolAz = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoRiseSolAlt = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoSetTime = cell(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoSetSolAz = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoSetSolAlt = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
diffuseRiseTime = cell(size(xyPoints,1),size(DOYs,2)); %character array
diffuseSetTime = cell(size(xyPoints,1),size(DOYs,2)); %character array
diffuseLength = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array

% Constants
numDays = length(DOYs);
numTimes = length(LocalTimes);
numPoints = size(xyPoints,1);

%% Calculate solar zenith and azimuth angles for a point
UTCtimes = LocalTimes - UTCoffset; %input is in local time, need to input times in function in UTC

for pnt = 1:size(xyPoints,1) % loop for each point in array
    
    % HemisphericalArea (Sky View Factor)
    SkyVF(pnt) = mean((90-illum_angle(pnt,:))./90); %from Yard paper
    % AzimuthOfLowestAltitude90to269 ?
    % AzimuthOfLowestAltitude270to89 ?
    % CanyonAzimuth ?
    % CanyonHemiAngle ?
    
    for jday = 1:numDays % loop for days of year to test

        DOY = DOYs(jday); %day of year for iteration
        % function to calculate zenith and azimuth angles from Margilus. Zenith angle are measured from the VERTICAL
        [solar_zenith,solar_azimuth,sunrise,sunset,solar_decl,hour_angle]...
            =solar_geometry(DOY,UTCtimes,UTCoffset,lat_deg(pnt),lon_deg(pnt)); %#ok<ASGLU> %supress warning

        %% Solor Geometry/ Topography/ Statistics from Yard Model
        % ----------- Solar Geometry -----------
        
        % SunRise Time
        sunriseTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(sunrise),'HH:MM:SS')); % string
        % SunRise Azimuth
        sunriseAz(pnt,jday) = interp1(LocalTimes,solar_azimuth,sunrise); %degrees
        % SunRise Zenith
        sunriseZe(pnt,jday) = interp1(LocalTimes,solar_zenith,sunrise); %degrees from vertical
        % TotalDayLength
        TotalDayLength(pnt,jday) = (sunset-sunrise)*60; % minutes
        % SunSet Time
        sunsetTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(sunset),'HH:MM:SS')); % string
        % SunSet Azimuth
        sunsetAz(pnt,jday) = interp1(LocalTimes,solar_azimuth,sunset); %degrees
        % SunSet Zenith
        sunsetZe(pnt,jday) = interp1(LocalTimes,solar_zenith,sunset); %degrees from vertical
        % DirectBeamDayLength
        solar_altitude = 90-solar_zenith; % degrees measured from the HORIZONTAL
        incidence(jday,:,pnt) = interp1(azimuthDeg,illum_angle(pnt,:),solar_azimuth) < solar_altitude; % logical
        directBeamDayLength(pnt,jday) = sum(diff(LocalTimes(incidence(jday,:,pnt))))*60; %minutes
        
        % ----------- Topography ---------------
        if directBeamDayLength(pnt,jday) ~= 0 %check to see if point gets no direct been
                            % TopoRise Time      
            TopoRiseT = hours(minutes(find(incidence(jday,:,pnt),1,'first'))); %finds the first 1 in the incidence array and calculates time duration in decimal hours
            TopoRiseTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(TopoRiseT),'HH:MM:SS')); % string
            % TopoRise solazimuth
            TopoRiseSolAz(pnt,jday) = interp1(LocalTimes,solar_azimuth,TopoRiseT); % degree
            % TopoRise solAltitude
            TopoRiseSolAlt(pnt,jday) = interp1(LocalTimes,solar_altitude,TopoRiseT); % degrees
            % TopoSet Time
            TopoSetT = hours(minutes(find(incidence(jday,:,pnt),1,'last'))); %finds the last 1 in the incidence array and calculates time duration in decimal hours
            TopoSetTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(TopoSetT),'HH:MM:SS')); % string
            % TopoSet solazimuth
            TopoSetSolAz(pnt,jday) = interp1(LocalTimes,solar_azimuth,TopoSetT); %degrees
            % TopoSet solAltitude
            TopoSetSolAlt(pnt,jday) = interp1(LocalTimes,solar_altitude,TopoSetT); % degrees from horizontal
        
        else %if no direct beam set values to zero
            % TopoRise Time      
            TopoRiseTime(pnt,jday) = {'NULL'};
            % TopoRise solazimuth
            TopoRiseSolAz(pnt,jday) = 0;
            % TopoRise solAltitude
            TopoRiseSolAlt(pnt,jday) = 0;
            % TopoSet Time
            TopoSetTime(pnt,jday) = {'NULL'};
            % TopoSet solazimuth
            TopoSetSolAz(pnt,jday) = 0;
            % TopoSet solAltitude
            TopoSetSolAlt(pnt,jday) = 0;
        end
        
        % ----------- Diffuse -----------
        % Diffuse Sunrise Time 
        diffuseRise = interp1(solar_zenith(1:floor(numTimes/2)),LocalTimes(1:floor(numTimes/2)),sunriseZe(pnt,jday)+6); %split arrays in half because zenith values repeat
        diffuseRiseTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(diffuseRise),'HH:MM:SS')); % string
        % Diffuse Sunset Time
        diffuseSet = interp1(solar_zenith(floor(numTimes/2):end),LocalTimes(floor(numTimes/2):end),sunsetZe(pnt,jday)+6); %split arrays in half because zenith values repeat
        diffuseSetTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(diffuseSet),'HH:MM:SS')); % string
        % Diffuse Length (Time)
        diffuseLength(pnt,jday) = (diffuseSet-diffuseRise)*60; % minutes

        % ----------- Insolation -----------
        % InsolationDirect ?
        % InsolationInDifuse ?
        % InsolationTotal ?       
        
    end
end

% reshape incidence to be 2d array where
% rows == datetimes and columns == xy points
% incidence = reshape(incidence,[numDays*numTimes,numPoints]);

%%
% Celestial variables:
for j = 4%:length(DOYs) %column
    for i = 1:10:191 %row
        disp({[SkyVF(i)];[];[];[];[];...
        char(sunriseTime(i,j));sunriseAz(i,j);TotalDayLength(i,j);char(sunsetTime(i,j));sunsetAz(i,j);...
        directBeamDayLength(i,j);char(TopoRiseTime(i,j));[];TopoRiseSolAz(i,j);...
        TopoRiseSolAlt(i,j);char(TopoSetTime(i,j));[];TopoSetSolAz(i,j);TopoSetSolAlt(i,j);...
        char(diffuseRiseTime(i,j));char(diffuseSetTime(i,j));diffuseLength(i,j);...
        [];[];[]})
        disp({[];[]})
    end
end
