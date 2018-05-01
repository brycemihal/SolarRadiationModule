function [SkyVF] = SVF(illum_angle) 
%% Function to calculate the sky view factor (SVF)
% Bryce Mihalevich
% 4/3/18

% Inputs:


% Outputs:

for i = 1:size(illum_angle,1)
    %% Change illumination angles to be from horizontal
    illum_angle_horiz = 90-illum_angle(i,:);
    
    %% Moore et. al. 2014 Method
    % from Johnson and Watson, 1984
    % just an average of illumination angles, where angles are measured from
    % the horizontal
    SkyVF(i).moore = mean(cosd(illum_angle_horiz));

    %% Yard et. al. 2005 Method
    % from Anton 1984; illumination angle measured from horizontal
    SkyVF(i).yard = mean((90-illum_angle_horiz)./90);

    %% Margulis 2017 MOD-WET toolbox
    % SVF from Dozier and Marks, 1987
    % Find average horizon angle (across all azimuth)
    mean_horizon_angle=mean(illum_angle_horiz);
    SkyVF(i).margulis=cosd(mean_horizon_angle).^2;

    %% from [http://dx.doi.org/10.1063/1.4921386]
    % Liu-Jordan method
    SkyVF(i).liujordan = mean((1+cosd(illum_angle_horiz))/2);
    % Tian et al
    SkyVF(i).tian = mean((180-illum_angle_horiz)/180);

end

%%
hold on
plot(cell2mat(({SkyVF.moore})))
plot(cell2mat(({SkyVF.yard})))
plot(cell2mat(({SkyVF.margulis})))
plot(cell2mat(({SkyVF.liujordan})))
plot(cell2mat(({SkyVF.tian})))
legend('Moore','Yard','Margulis','Liu-Jordan','Tian')
%% References
% Margulis
% [http://dx.doi.org/10.1063/1.4921386]

% Moore, R. D., Leach, J. A., & Knudson, J. M. (2014). Geometric calculation of view factors for stream surface radiation modelling in the presence of riparian forest. Hydrological Processes, 28(6), 2975?2986. https://doi.org/10.1002/hyp.9848
% Yard, M. D., Bennett, G. E., Mietz, S. N., Coggins, L. G., Stevens, L. E., Hueftle, S., & Blinn, D. W. (2005). Influence of topographic complexity on solar insolation estimates for the Colorado River, Grand Canyon, AZ. Ecological Modelling, 183(2?3), 157?172. https://doi.org/10.1016/j.ecolmodel.2004.07.027