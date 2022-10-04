function [Dynamics] = SyncLidarData(input, curFileInfo, Dynamics, windfield)
%This function extends the dynamics time series to account for the
%measurement distance. It alligns the measurement times defined in the
%lidar configuration with the dynamics time series and extracts required
%time steps.

%OUTPUTS
%-Dynamics Struct, containing the extended dynamics timeseries of the
%required measurement times according to the configured lidar measurement
%timesteps
%
% V.Pettas/M.Gräfe
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021

%%
fNames = fieldnames(Dynamics.sim.channels);

beams = length(input.PatternY{cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Pat')))),1});
% Get sampling rate and total simulation time
Dynamics.sim.dt = Dynamics.sim.channels.(fNames{10})(2,1)-Dynamics.sim.channels.(fNames{10})(1,1);


%% Extend dynamics time series to account for measurement distance
ref_plane_dist = input.ref_plane_dist;
URef = windfield.URef;
dt = Dynamics.sim.dt;

% find number of samples which need to added to account for measurement Distance
num = floor((ref_plane_dist/ URef)/dt);
Dynamics.num_chunk = num;

for ch = 1:length((fNames))
    if num > 2 % only do the extension if more than two samples are added
        chunk = Dynamics.sim.channels.(fNames{ch})(num:2*num-1,1); % find chunk of channel
        Dynamics.sim.channels.(fNames{ch}) = [chunk;Dynamics.sim.channels.(fNames{ch})]; % add in the beginning
        % Smooth transition by replacing with moving  vaverage
        num_av = 10;% find number of points to be averaged on each side of transition based on time of 5s. Make this a parameter?
        if num<num_av % if chunk is smaller than intended number replace
            num_av = num-1;
        end
        Dynamics.sim.channels.(fNames{ch})((num-num_av:num+num_av),1) = movmean(Dynamics.sim.channels.(fNames{ch})(num-num_av:num+num_av,1),(num_av/2));
    end
    Dynamics.sim.channels.(fNames{ch}) = repmat(Dynamics.sim.channels.(fNames{ch})',beams,1); % create matrix with one row per beam
    
end

% resample time to account for extended dynamics
for i = 1:beams
    Dynamics.sim.channels.(fNames{10})(i,:) = (0:(length(Dynamics.sim.channels.(fNames{10})(i,:))-1))*dt;
end

%overwrite DOF with zero depending on input.DOF.tilt, input.DOF.trans
channelnames = input.channels.channelnames;
if input.DOF.rot_trans == false
    Dynamics.sim.channels.(channelnames{1})(:,:) = 0;
    Dynamics.sim.channels.(channelnames{2})(:,:) = 0;
    Dynamics.sim.channels.(channelnames{3})(:,:) = 0;
    
    Dynamics.sim.channels.(channelnames{7})(:,:) = 0;
    Dynamics.sim.channels.(channelnames{8})(:,:) = 0;
    Dynamics.sim.channels.(channelnames{9})(:,:) = 0;
end

if input.DOF.vel == false
    Dynamics.sim.channels.(channelnames{4})(:,:) = 0;
    Dynamics.sim.channels.(channelnames{5})(:,:) = 0;
    Dynamics.sim.channels.(channelnames{6})(:,:) = 0;
end

%% get total sim length

fns = fieldnames(Dynamics.sim.channels);
Dynamics.sim.sim_time = length(Dynamics.sim.channels.(fns{1}));

%% create x_L,y_L,z_L matrix
timeStep_Measurements = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Tm'))));  % [s] Time required for a single point measurement
timestep_pat_vec = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Tp'))));     % [s] Time required for a full scan


fullTime         = (windfield.dt*windfield.grid.nt)-windfield.dt;
num_measurements = length(floor((1:(size(windfield.u,2)*(timestep_pat_vec/fullTime)):windfield.grid.nt))); % find number of measurements fitting into simulated timeseries

% Create matrices of lidar pattern for each scan
y_L = input.PatternY{cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Pat')))),1}; % lidar pattern coordinates lateral (m)
z_L = input.PatternZ{cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Pat')))),1}; % lidar pattern coordinates vertical (m)
x_L(1,1:length(y_L)) = input.ref_plane_dist;

y_L = y_L'-input.Pos_LiDAR(1);
z_L = z_L'-input.Pos_LiDAR(2);
x_L = x_L';

x_L = repmat(x_L,[1 num_measurements]);
y_L = repmat(y_L,[1 num_measurements]);
z_L = repmat(z_L,[1 num_measurements]);

Dynamics.x_L = x_L;
Dynamics.y_L = y_L;
Dynamics.z_L = z_L;


%find timestamps where lidar measurements take place and extract only these
for ch = 1:length((fNames))-1 % all channels besides time
    for i = 1:beams
        time_Measurement=((i-1)*timeStep_Measurements : timestep_pat_vec : fullTime);
        % find correct value in dynamics time series by 'nearest'
        % interpolation based on measurement time per beam
        channel_temp = interp1(Dynamics.sim.channels.(channelnames{10})(i,:),Dynamics.sim.channels.(fNames{ch})(i,:),time_Measurement,'nearest');
        if length(channel_temp)<num_measurements % extend with NaN if necessary
            channel_temp(1,end+1:num_measurements) = NaN;
        end
        channel_temp = channel_temp(1,1:num_measurements); % cut to number of measurements if necessary
        Dynamics.sim.channels_time.(fNames{ch})(i,:) = channel_temp; % store to struct
        clear channel_temp
    end
end
end
