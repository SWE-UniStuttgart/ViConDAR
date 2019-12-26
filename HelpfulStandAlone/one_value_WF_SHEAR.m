%% HEADER

% This code gets an .hh file and creates a .mat windfield (Turbsim style) with
% constant wind in every slice based on the .hh file. Additionally also shear can be added
% to the wind file. V and W components are considered 0. 

% Francisco Costa
% © SWE
%obtain Windfield data 
%Inputs
%Outputs
%Dependencies

clc;
close all;
clearvars

%--------------------------------------------------------------------------%

%% Input parameter definition:

save_WF_directory  = '../Output/'; % directory to save the output .mat file
FileName           = 'testWindfield_noshear'; % file name to be saved the variable will always be windfield

hh_input  = '../WindFiles/DTU10MW_DLC12_ClassA_1_04_1.hh';  %input hh file to readin speeds and time

%requested grid
ny  = 41; %define amount of point in y direction of grid
nz  = 41; %define amount of point in z direction of grid
dy  = 3.25; %define spacing in y direction of grid
dz  = 3.25; %%define spacing in z direction of grid
gridLimY = 65; %grid limit in y direction assuming [-gridLimY : gridLimY]
gridLimZ = 65; %grid limit in z direction assuming [-gridLimY : gridLimY] 
hub_h    = 119;

%Shear
shear_alpha = 0.0; %Apply vertical shear if zero there will be uniform wind in each time slice

%% Read in hh file

HH_data=importdata(hh_input);

%% Start populating windfield with known values
windfield.ny = ny;
windfield.nz = nz;
windfield.dt = HH_data.data(2,1)-HH_data.data(1,1);
windfield.URef = round (mean(HH_data.data(:,2)));
windfield.T_offset = 0;
windfield.grid.nt = length(HH_data.data(:,2));
windfield.grid.dy = dy;
windfield.grid.ny = ny;
windfield.grid.nz = nz;
windfield.grid.dt = windfield.dt;
windfield.grid.t  = HH_data.data(:,1);
windfield.grid.dz = dz;
windfield.grid.z = -gridLimZ:dz:gridLimZ;
windfield.grid.y = -gridLimY:dz:gridLimY;
windfield.grid.Y = repmat(windfield.grid.y,nz,1) ;
windfield.grid.Z = repmat(windfield.grid.z',1,ny) ;
windfield.v      = zeros (ny,length(HH_data.data(:,2)),nz); % v+w are considered zero in this implementation
windfield.w      = zeros (ny,length(HH_data.data(:,2)),nz); % v+w are considered zero in this implementation
windfield.u      = zeros (ny,length(HH_data.data(:,2)),nz); % preallocation

%Create vertical shear for each slice and put in the correct slices 
for i = 1:windfield.grid.nt
    u_vec(i,:) = HH_data.data(i,2).*( (windfield.grid.z+hub_h)/hub_h).^shear_alpha; %#ok<SAGROW>
    windfield.u (:,i,:) = repmat(u_vec(i,:),ny,1);
end
save([save_WF_directory FileName],'windfield')


