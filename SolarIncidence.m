function[incidence,dateTimes,SkyVF] = SolarIncidence(xyPoints,xyCoord,azimuthDeg,illum_angle,DOYs,LocalTimes,UTCoffset)
% Bryce Mihalevich
% Last Modified: 9/6/18
%
% Description: function to determine solar incidence and sky view factor.
% Several other output variables are calculated but are not included in the
% output of the function at this time. If certain variables are needed
% add them to the brackets above to include in output function. Reminder, 
% variables need to change where it function is called too if this is done. 
%
% Input Variables: 
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
% Potential Output Variables: add variables to function brackets for output
% ------Shade variables:
%       incidence = binary array indicating when a point recieves direct
%       beam (DB) radiation (1's) and when it does NOT recieve DB (0's).
%       Row in array equals datetimes (LocalTimes for each DOYs)
%       Columns in array equal number of xy points
%       Future: 3rd dim for cross section width locations with each xyPoint
%       dateTimes = datetime array with LocalTimes for each DOY
% 
% ------Constant with time variables:
%       SkyVF = Sky view factor or hemispherical area
%       AzimuthOfLowestAltitude90to269 = unknown/not calculated
%       AzimuthOfLowestAltitude270to89 = unknown/not calculated
%       CanyonAzimuth = unknown/not calculated
%       CanyonHemiAngle = unknown/not calculated
%
% ------Solar Position variables:
%       sunriseTime = calculated time of sunrise specific to a day and point
%       sunriseAz = azimuth angle at sunrise specific to a day and point
%       sunriseZe = zenith angle at sunrise specific to a day and point
%       sunsetTime = calculated time of sunset specific to a day and point
%       sunsetAz = azimuth angle at sunset specific to a day and point
%       sunsetZe = zenith angle at sunset specific to a day and point
%       TotalDayLength = sunrise minus sunset in minutes
%       
% ------Topo variables:
%       TopoRiseTime = time when sun first hist point of interest  
%       TopoRiseSolAz = azimuth angle when sun first hits point of interst
%       TopoRiseSolAlt = 90-zenith angle when sun first hits point of interest
%       TopoRiseIncidence = unknown/not calculated
%       TopoSetTime = time when sun first hist point of interest  
%       TopoSetSolAz = azimuth angle when sun first hits point of interst
%       TopoSetSolAlt = 90-zenith angle when sun first hits point of interest
%       directBeamDayLength = TopoSetTime minut TopoRiseTime in minutes
%       TopoSetIncidence = unknown/not calculated
%       
% ------Radiation variabels:
%       diffuseRiseTime = time of morning civil twilight, defined as 6 degrees
%       (from the vertical) below the horizon (i.e.: ~96 degrees)
%       diffuseSetTime = time of evening civil twilight, defined as 6 degrees
%       (from the vertical) below the horizon  (i.e.: ~96 degrees)
%       difuseLength = diffusesetTime minus diffuseRiseTime in minutes
%       InsolationDirect = unknown/not calculated
%       InsolationDiffuse = unknown/not calculated
%       InsolationTotal = unknown/not calculated
%
% For example code see: Example_Script_to_use_functions.m in Examples folder. 
%
disp('running SolarIncidence.m');
%% Determine coordinate system of points and convert xyCoords if needed
% if xy coordinates are in UTM then convert
if strcmp(xyCoord,'UTM') %then convert to UTM
        zone = 12;
        isnorth = 1;
        [lat,lon] = utmups_inv(xyPoints(:,1),xyPoints(:,2),zone,isnorth); %UTM to DD conversion
        xyPoints = [lon,lat];
end

% separate x and y values
lon_deg = xyPoints(:,1);
lat_deg = xyPoints(:,2);

%% preallocate memory to variables for speed
% Arrays
SkyVF = zeros(size(xyPoints,1),1); %floating point array
southHemiMinAz = zeros(size(xyPoints,1),1); %floating point array
northHemiMinAz = zeros(size(xyPoints,1),1); %floating point array
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
TopoRiseSolZe = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoRiseSolAz = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoRiseSolAlt = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoSetTime = cell(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoSetSolZe = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoSetSolAz = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
TopoSetSolAlt = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array
diffuseRiseTime = cell(size(xyPoints,1),size(DOYs,2)); %character array
diffuseSetTime = cell(size(xyPoints,1),size(DOYs,2)); %character array
diffuseLength = zeros(size(xyPoints,1),size(DOYs,2)); %floating point array

% Constants
numDays = length(DOYs);
numTimes = length(LocalTimes);
numPoints = size(xyPoints,1);
southernHemAngles = 90:269;
northernHemAngles = [0:89,270:359];
%% Calculate solar zenith and azimuth angles for a point
UTCtimes = LocalTimes - UTCoffset; %input is in local time, need to input times in function in UTC

for pnt = 1:size(xyPoints,1) % loop for each point in array
    
    % HemisphericalArea (Sky View Factor)
    SkyVF(pnt) = mean((90-illum_angle(pnt,:))./90); %from Yard paper - consider changing to other svf methods
%     % AzimuthOfLowestAltitude90to269 = lowest altitude angle in southern hemisphere 
%     [~,azimuth] = min(illum_angle(pnt,90:269));
%     southHemiMinAz(pnt) = southernHemAngles(azimuth);
%     % AzimuthOfLowestAltitude270to89 = lowest altitude angle in northern hemisphere
%     [~,azimuth] = min(illum_angle(pnt,[1:89,270:end]));
%     northHemiMinAz(pnt) = northernHemAngles(azimuth);
%     % CanyonAzimuth = quick way to estimate direction of the river
%     % CanyonHemiAngle =  used to get  HemisphericalArea (Sky View Factor) for diffuse insolation

    for jday = 1:numDays % loop for days of year to test

                DOY = DOYs(jday); %day of year for iteration
        % function to calculate zenith and azimuth angles from Margilus.
        % Zenith angle are measured from the VERTICAL
        % arrays are the length of num times
        [solar_zenith,solar_azimuth,sunrise,sunset,solar_decl,hour_angle]...
            =solar_geometry(DOY,UTCtimes,UTCoffset,lat_deg(pnt),lon_deg(pnt)); %#ok<ASGLU> %supress warning

        %% Solor Geometry/ Topography/ Statistics from Yard Model
        % ----------- Solar Geometry -----------
        
%         % SunRise Time
%         sunriseTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(sunrise),'HH:MM:SS')); % string
%         % SunRise Azimuth
%         sunriseAz(pnt,jday) = interp1(LocalTimes,solar_azimuth,sunrise); %degrees
%         % SunRise Zenith
%         sunriseZe(pnt,jday) = interp1(LocalTimes,solar_zenith,sunrise); %degrees from vertical
%         % TotalDayLength
%         TotalDayLength(pnt,jday) = (sunset-sunrise)*60; % minutes
%         % SunSet Time
%         sunsetTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(sunset),'HH:MM:SS')); % string
%         % SunSet Azimuth
%         sunsetAz(pnt,jday) = interp1(LocalTimes,solar_azimuth,sunset); %degrees
%         % SunSet Zenith
%         sunsetZe(pnt,jday) = interp1(LocalTimes,solar_zenith,sunset); %degrees from vertical
        % DirectBeamDayLength
        solar_altitude = 90-solar_zenith; % degrees measured from the HORIZONTAL
        incidence(jday,:,pnt) = interp1(azimuthDeg,illum_angle(pnt,:),solar_azimuth) < solar_altitude; % logical
        directBeamDayLength(pnt,jday) = sum(diff(LocalTimes(incidence(jday,:,pnt))))*60; %minutes
        
        % ----------- Topography ---------------
%         if directBeamDayLength(pnt,jday) ~= 0 %check to see if point gets no direct been
%             % TopoRise Time      
%             TopoRiseT = hours(minutes(find(incidence(jday,:,pnt),1,'first'))); %finds the first 1 in the incidence array and calculates time duration in decimal hours
%             TopoRiseTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(TopoRiseT),'HH:MM:SS')); % string
%             % TopoRise solar Zenith (incidence)
%             TopoRiseSolZe(pnt,jday) = interp1(LocalTimes,solar_zenith,TopoRiseT); % degree
%             % TopoRise solar Azimuth
%             TopoRiseSolAz(pnt,jday) = interp1(LocalTimes,solar_azimuth,TopoRiseT); % degree
%             % TopoRise solar Altitude
%             TopoRiseSolAlt(pnt,jday) = interp1(LocalTimes,solar_altitude,TopoRiseT); % degrees
%             % TopoSet Time
%             TopoSetT = hours(minutes(find(incidence(jday,:,pnt),1,'last'))); %finds the last 1 in the incidence array and calculates time duration in decimal hours
%             TopoSetTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(TopoSetT),'HH:MM:SS')); % string
%             % TopoSet solar Zenith (incidence)
%             TopoSetSolZe(pnt,jday) = interp1(LocalTimes,solar_zenith,TopoSetT); %degrees
%             % TopoSet solar Azimuth
%             TopoSetSolAz(pnt,jday) = interp1(LocalTimes,solar_azimuth,TopoSetT); %degrees
%             % TopoSet solar Altitude
%             TopoSetSolAlt(pnt,jday) = interp1(LocalTimes,solar_altitude,TopoSetT); % degrees from horizontal
%         
%         else %if no direct beam set values to zero
%             % TopoRise Time      
%             TopoRiseTime(pnt,jday) = {'NULL'};
%             % TopoRise solazimuth
%             TopoRiseSolAz(pnt,jday) = 0;
%             % TopoRise solAltitude
%             TopoRiseSolAlt(pnt,jday) = 0;
%             % TopoSet Time
%             TopoSetTime(pnt,jday) = {'NULL'};
%             % TopoSet solazimuth
%             TopoSetSolAz(pnt,jday) = 0;
%             % TopoSet solAltitude
%             TopoSetSolAlt(pnt,jday) = 0;
%         end
%         
%         % ----------- Diffuse -----------
%         % Diffuse Sunrise Time 
%         diffuseRise = interp1(solar_zenith(1:floor(numTimes/2)),LocalTimes(1:floor(numTimes/2)),sunriseZe(pnt,jday)+6); %split arrays in half because zenith values repeat
%         diffuseRiseTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(diffuseRise),'HH:MM:SS')); % string
%         % Diffuse Sunset Time
%         diffuseSet = interp1(solar_zenith(floor(numTimes/2):end),LocalTimes(floor(numTimes/2):end),sunsetZe(pnt,jday)+6); %split arrays in half because zenith values repeat
%         diffuseSetTime(pnt,jday) = strcat(datestr(DOY,'dd-mmm'),{' '},datestr(hours(diffuseSet),'HH:MM:SS')); % string
%         % Diffuse Length (Time)
%         diffuseLength(pnt,jday) = (diffuseSet-diffuseRise)*60; % minutes
% 
%         % ----------- Insolation -----------
%         % not sure how to calculate these yet
%         % InsolationDirect ?
%         % InsolationInDifuse ?
%         % InsolationTotal ?       
        
    end
    % Percent Complete
    if mod(pnt,round(size(xyPoints,1)/50))==0 %write out percent progress to command window
        percentComplete = round(pnt/size(xyPoints,1)*100);
        fprintf('~%d%% complete\n',percentComplete);
    end
end

% reshape incidence to be 2d array where
% rows == datetimes and columns == xy points
incidence = reshape(incidence,[numDays*numTimes,numPoints]);

%% Build datetime array
%reshape DOYs to be number of LocalTimes for each DOY
dateArray = repmat(DOYs,size(LocalTimes,2),1); %repeat matrix for number of values in LocalTimes
dateArray = dateArray(:); %stack columns
timeArray = repmat(LocalTimes',size(DOYs,2),1); %repeat matrix for number of DOYs

datetimeString = strcat(datestr(datestr(dateArray,'dd-mmm'),'yyyy-mm-dd'),{' '},datestr(hours(timeArray),'HH:MM:SS')); %concatenate the date and the time
dateTimes = datetime(datetimeString,'InputFormat','yyyy-MM-dd HH:mm:SS'); %auto appends the current year
