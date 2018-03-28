function [zenith_angle_deg,azimuth_angle_deg,sunrise,sunset,solar_decl,hour_angle]...
    =solar_geometry(DOY,UTC,time_zone_shift,lat_deg,lon_deg)
% -------------------------------------------------------------------------
% [zenith_angle_deg,azimuth_angle_deg,sunrise,sunset,solar_decl,hour_angle]...
%     =solar_geometry(DOY,UTC,time_zone_shift,lat_deg,lon_deg)
% 
% Description:
% Function to compute solar zenith and azimuth angles, hour of sunrise and
% sunset, and the solar declination and hour angles
%
% Written by Steve Margulis, 6/20/2011
%
% Inputs:
% DOY: Day of Year (integer): Jan. 1 = 1, Dec. 31 = 365
% UTC: Universal Time (Greenwich Mean Time) in decimal hours: 0.0=midnight, 23.0=11pm
% time_zone_shift (hours): Time shift relative to UTC (hours) (i.e. California, -8)
% lat_deg: Latitude (in degrees)
% lon_deg: Longitude (in degrees), negative to West of UTC line
% 
% Outputs: 
% zenith angle_deg: Solar zenith angle in degrees
% azimuth_angle_deg: Solar azimuth angle in degrees
% sunrise: Local hour of sunrise
% sunset: Local hour of sunset
% solar_decl: Solar declination angle in radians
% hour_angle: Hour angle in radians
%
% NOTE: Can run for multiple hours (i.e. vector inputs), but for a single
% day (DOY) and lat/lon
%
% Example:
% DOY= 115;% April 25
% UTC= 16;
% time_zone_shift= -8;
% lat_deg= 34.05; lon_deg= -118.25; % Los Angeles
% 
% [zenith_angle_deg,azimuth_angle_deg,sunrise,sunset]=...
%     solar_geometry(DOY,UTC,time_zone_shift,lat_deg,lon_deg);
% 
% zenith_angle_deg =
%    56.2121
% azimuth_angle_deg =
%    97.3170
% sunrise =
%     5.2555
% sunset =
%    18.4450
% solar_decl =
%     0.2256
% hour_angle =
%    -1.0080

%%% Disclaimer:
%%% This program and all related codes that are part of the "MOdular 
%%% Distributed Watershed Educational Toolbox" (hereafter "MOD-WET" or 
%%% "software") is designed for instructional and educational use. 
%%% Commercial use is prohibited. The software is provided 'as is' without 
%%% warranty of any kind, either express or implied. MOD-WET could include 
%%% technical or other mistakes, inaccuracies or typographical errors. The 
%%% use of the software is done at your own discretion and risk and with 
%%% agreement that you will be solely responsible for any damage and that 
%%% the authors and their affiliate institutions accept no responsibility 
%%% for errors or omissions in the software or documentation. In no event 
%%% shall the authors or their affiliate institutions be liable to you or 
%%% any third parties for any special, indirect or consequential damages of 
%%% any kind, or any damages whatsoever.
%%%
%%% Any bugs that are found can be reported to: Steve Margulis
%%% (margulis@seas.ucla.edu) where we will make every effort to fix them 
%%% for future releases.

latrad = lat_deg * pi / 180;

time_local=UTC+time_zone_shift;
if time_local<0
    time_local=24+time_local;
    DOY=DOY-1;
end

day_angle = 2. * pi * (DOY - 1.) / 365.; % day angle in radians
% Declination angle in radians
solar_decl =  0.006918 - 0.399912 .* cos(day_angle) +                     ...
      0.070257 .* sin(day_angle) - 0.006758 .* cos(2.* day_angle) +...
      0.000907 .* sin(2. * day_angle) -                           ...
      0.002697 .* cos(3. * day_angle) + 0.00148 .* sin(3. * day_angle);

LSTM=15*time_zone_shift; % local standard time meridian (deg.) -- Time zones defined by 15 deg. longitude increments around the globe
B=360/365*(DOY-81); % degrees
% Equation of Time
EofT_min2=9.87*sin(2*B*pi/180) - 7.53*cos(B*pi/180) -1.5*sin(B*pi/180);
TC=4*(lon_deg-LSTM)+EofT_min2;  % time correction (minutes)
LST=time_local +TC/60;  % local solar time (hours)
hour_angle=15*(LST-12.)*pi/180.;  % radians

%Compute the solar zenith angle (degrees)
zenith_angle = acos(sin(latrad) * sin(solar_decl) + cos(latrad) * cos(solar_decl) .* cos(hour_angle));

% Compute solar azimuth angle (degrees)
azimuth_angle= acos((sin(solar_decl).*cos(latrad)-cos(solar_decl).*sin(latrad).*cos(hour_angle))./sin(zenith_angle));
azimuth_angle(LST>12)=2*pi-azimuth_angle(LST>12);

zenith_angle_deg=zenith_angle*180/pi;
azimuth_angle_deg=azimuth_angle*180/pi;

% Sunrise/sunset at local time
sunrise=12-180/15/pi*acos(-sin(latrad)*sin(solar_decl)./cos(latrad)./cos(solar_decl))-TC/60.;
sunset=12+180/15/pi*acos(-sin(latrad)*sin(solar_decl)./cos(latrad)./cos(solar_decl))-TC/60.;

return