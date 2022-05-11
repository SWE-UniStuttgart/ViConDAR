%% Header
%
% This Script can be used to read Output from aeroelastic Simulation (eg.
% openFAST) and store dyanmics into dynamics input struct for ViConDAR.
% This is an specific implementation for the SWE openFAST post porcessing
% work flow. Needs tobe adapted to the users needs. 
% M.Gräfe
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021


clc
clear all
% set inputs
timetocut = 600; % time that is cut off in the begining of timesereries to remove transition phase
FAST_directory = 'FAST_Results\'; % Directory of FAST outputs. All files in this directory will be processed
input.OriginalDynamics_dir = 'Original_Dynamics'; % Directory to store input dyanmics for use in ViConDAR
filenames = dir(fullfile(FAST_directory,'*.mat')); % List of files available in directory
[n_i,~] = size(filenames);

%% input channels from FAST

% Dynamic channel names from FAST
input.channels.pitch='PtfmPitch';% Pitch rotational displacement, Format: Outout.sim.channels.channelnames
input.channels.roll='PtfmRoll';% Roll rotational displacement, Format: Outout.sim.channels.channelnames
input.channels.yaw='PtfmYaw';% Yaw rotational displacement, Format: Outout.sim.channels.channelnames

input.channels.x_vel_L='NcIMUTVxs'; % translational velociy of Lidar in Lidar CS
input.channels.y_vel_L='NcIMUTVys'; % translational velociy of Lidar in Lidar CS
input.channels.z_vel_L='NcIMUTVzs'; % translational velociy of Lidar in Lidar CS

input.channels.x_trans_L='TwrTpTDxi';% heave in Lidar CS, Format: Outout.sim.channels.channelnames
input.channels.y_trans_L='TwrTpTDyi';% surge in Lidar CS, Format: Outout.sim.channels.channelnames
input.channels.z_trans_L='TwrTpTDzi';% sway in Lidar CS, Format: Outout.sim.channels.channelnames

%%Output Channelnames for vicondar
input.outchannels.pitch='Pitch';% Pitch rotational displacement, Format: Outout.sim.channels.channelnames
input.outchannels.roll='Roll';% Roll rotational displacement, Format: Outout.sim.channels.channelnames
input.outchannels.yaw='Yaw';% Yaw rotational displacement, Format: Outout.sim.channels.channelnames

input.outchannels.x_vel_L='x_vel'; % translational velociy of Lidar in Lidar CS
input.outchannels.y_vel_L='y_vel'; % translational velociy of Lidar in Lidar CS
input.outchannels.z_vel_L='z_vel'; % translational velociy of Lidar in Lidar CS

input.outchannels.x_trans_L='x_trans_L';% heave in Lidar CS, Format: Outout.sim.channels.channelnames
input.outchannels.y_trans_L='y_trans_L';% surge in Lidar CS, Format: Outout.sim.channels.channelnames
input.outchannels.z_trans_L='z_trans_L';% sway in Lidar CS, Format: Outout.sim.channels.channelnames


for i_file=1:n_i
  
    % call functions
    Output = struct();
    cd (FAST_directory);
    load(filenames(i_file).name); % load input file from FAST, with variable TimeResults
    Output.sim.TimeResults = TimeResults; 
 
    %% sort channels into struct with output fieldnames
    fNames=fieldnames(Output.sim.TimeResults);
    finNames=fieldnames(input.channels);
    fOutNames=fieldnames(input.outchannels);
    dt=TimeResults.Time{1, 1}(2,1)-TimeResults.Time{1, 1}(1,1);
    num_cut=timetocut/dt;
    
    for i=3:length(fNames)
        for f=1:length(fOutNames)
            if strcmp(fNames{i},input.channels.(finNames{f}))
                x=Output.sim.TimeResults.(fNames{i}).Data;
                % cut channel to desired length
                x{1, 1}=x{1,1}(num_cut+1:end,1);
                
                Output.sim.channels.(input.outchannels.(fOutNames{f}))=x{1,1};
            end
        end
    end
    % add sim time to struct
    Output.sim.channels.Time=Output.sim.TimeResults.Time{1,1}(num_cut+1:end,1);
    
    %% Correct offset for translational displacement channels
    input.PositionTrans_x=0; % 0 for x_i
    input.PositionTrans_y=0; % 0 for y_i
    input.PositionTrans_z=0; % 0 for z_i, in FAST this is relative to MSL !!!
    
    Output.sim.channels.x_trans_L=Output.sim.channels.x_trans_L-input.PositionTrans_x;
    Output.sim.channels.y_trans_L=Output.sim.channels.y_trans_L-input.PositionTrans_y;
    Output.sim.channels.z_trans_L=Output.sim.channels.z_trans_L-input.PositionTrans_z;
    
    %% save Timeresults
    [pathstr,name,ext] = fileparts(filenames(i_file).name);
    Outname= strcat(name, '_dynamics');
    
    fastpath= input.OriginalDynamics_dir;
    cd (fastpath);
    new_name=strcat(Outname, '.mat');
    Dynamics=Output.sim.channels;
    save(new_name, 'Dynamics');
    
end



