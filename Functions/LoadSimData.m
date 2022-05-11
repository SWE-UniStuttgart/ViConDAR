function [Dynamics] = LoadSimData(directory,curFileInfo)

%  This function loads the timeseries input files 
% (see example format in TestFiles ('FLTEST_SD01_V10_TI09_dynamics.mat') the struct name should be 'Dynamics')
% and tranforms it to the interanl vicondar structure
%
%OUTPUT: 
% timeresult data stored in Output struct
% V.Pettas/M.Gräfe
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021
         

filename = strcat(curFileInfo.originalWF{1},'_dynamics');
PathToDYN = fullfile([directory,filename,'.mat']);

Dynamics_temp = load(PathToDYN); % Dynamics variable loaded
Dynamics = Dynamics_temp.Dynamics;

end

