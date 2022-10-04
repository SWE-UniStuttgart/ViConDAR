%% Header
%
% Lidar simulator. Create the output file including timeseries, errors and
% statistics from the lidar measurements.
%
% V.Pettas/F.Costa/M.Gräfe
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODOs:
% 1) Add the possiblility to rotate the field for yawed and inclined inflow
% 2) Check again the offset in combination with LOS
% 3) Handling lenght of lidar time series in cases first entry not existing  due to dynamics
% 4) consider to move Statistics, Resampling functions, noise adding, Shear calculation, REWS calculation out of getLidarOutput
% 5) consider option for shear calculation only for lidar measurement times to save time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------

function Output = getLidarOutput(input,curFileInfo)


%% IO/parameter definition

% Turbine parameters
rotor_radius = input.rotor_radius; % Radio of the Rotor [m]
Zh           = input.Zh;           % HubHeight [m]
Pos_LiDAR    = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Pos')))); %#ok<*FNDSB> % LiDAR position offsetfrom hub center(meters)==> [Y,Z] WE NEED TO FIX THE APPLICATION OF OFFSET IN LOS ONLY!!!!
% Lidar parameters
ref_plane_dist = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Fd'))));       % Reference Plane for LOS (distance[m])
distance_av_space = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'DAv'))));   % [m] values to use for imitating range gate averaging in the calcualtion of wind speeds from pulses meters ahead and afer the point
points_av_slice   = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'SlAv'))));  % how many point/slices you want to take in the averaging of distance_av_slice  Totalpoints = distance_av_slice/points_av_slice+1 IT HAS TO BE AN EXACT DIVISION FOR NOW!!!!
Y = input.PatternY{cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Pat')))),1}; % lidar pattern coordinates lateral (m)
Z = input.PatternZ{cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Pat')))),1}; % lidar pattern coordinates vertical (m)
timeStep_Measurements = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Tm'))));% [s] Time required for a single point measurement
timestep_pat_vec = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Tp'))));     % [s] Time required for a full scan

if timestep_pat_vec < timeStep_Measurements*length(Y) % throw error if timesteps dont match
    error ('Time step of pattern is smaller than the timestep of measures x number of points. Cannot continue')
end

%Processing REWS:
dist_REWS_nd = input.dist_REWS_nd; % Non dimentional span position for rotor effective wind speed calculation [define from 0 to 1 inclusive]
Wi           = input.Wi;  %Weight to be applied for rotor effective wind speed calculation
%Resampling lidar measurements, currently in Frequency domain
resampling_factor = input.resampling_factor; % Amount of desired resampling for outputs in Turbsim and PyConTurb used with input.flag_resampling
% Wind velocity vector components to be used
nComp             = input.nComp;        %#ok<NASGU> %1:u, 2:v+u 3:u+v+w % Number of components to process (U,V,W):
% Interpolation for time and space
type_interpolation     = input.type_interpolation; %#ok<NASGU> % (interp1) interpolation between slices [Time] line460 (check other options of interpm)
type_interpolation_2   = input.type_interpolation_2; % (interp2)  interpolation in grid points for values on the pattern points[Space]
%Noise magnitude to imitate uncertainty and noise in real measurements (in dB)
noise_U  = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Ns'))));% magnitude of noise to be applied in U time series (see help of awgn function)
noise_V  = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Ns')))); % magnitude of noise to be applied in V time series (see help of awgn function)
noise_W  = cell2mat(curFileInfo.values (find(strcmp((curFileInfo.variables{1, 1}),'Ns')))); % magnitude of noise to be applied in W time series (see help of awgn function)

%% Load windfield file
filenameOrWF = [input.OriginalWF_dir curFileInfo.originalWF{1} '.mat'];
load (filenameOrWF); %#ok<*LOAD> % winfield variable loaded

%% Get Dynamic Pattern from Aero Elastic Simulation

if input.flag_floating
    
    Pattern_n = length(Y);
    % Load Floater Dynamics Files for current wind file
    Dynamics = LoadSimData(input.OriginalDynamics_dir, curFileInfo);
    % Synchronize loaded floater dynamics with measurement times and extract rquired timestamps
    Dynamics = SyncLidarData(input, curFileInfo, Dynamics, windfield);
    % convert positions of focus points from lidar coordinate system to earth fixed coordinate system by applying rotation matrix based on dynamics
    Dynamics = Lidar2InertialDynamic(input, Dynamics);
    % extract Lidar position, translations, and translational velocities from dynamics input channels
    Dynamics = LidarTranslations(input, Dynamics);
    
    %%%%%%%%%%%% create variables for further LiDAR simulation
    x_I = cell(1,Pattern_n);
    y_I = cell(1,Pattern_n);
    z_I = cell(1,Pattern_n);
    
    for i = 1:Pattern_n
        x_I{1,i} = Dynamics.x_I(i,:);
        y_I{1,i} = Dynamics.y_I(i,:);
        z_I{1,i} = Dynamics.z_I(i,:);
        
        y_posL_I{1,i} = Dynamics.y_posL_I(i,:);
        z_posL_I{1,i} = Dynamics.z_posL_I(i,:);
        
        x_vel{1,i} = Dynamics.x_vel(i,:);
        y_vel{1,i} = Dynamics.y_vel(i,:);
        z_vel{1,i} = Dynamics.z_vel(i,:);
    end
    
    %assign Dynamics struct to output
    if input.flag_saveDynamicstoOutput
        Output.Dynamics = Dynamics;
    end
end

%% Obtain and create data
%extract components from the windfield variable

dt                  =  windfield.dt; %time step
gridtime            =  windfield.grid.nt;
gridny              =  windfield.grid.ny;
gridnz              =  windfield.grid.nz;
gridz               =  -windfield.grid.z;
gridy               =  windfield.grid.y;
Uref                =  windfield.URef; % Mean velocity of the windfied (m/s)
dz                  =  windfield.grid.dz;
dy                  =  windfield.grid.dy;
distanceSlices      =  Uref*dt; % Distance step  between consecutive slices(m)

% Manipulation of data before calculations:
compU = zeros(gridnz,gridtime,gridny);
compV = zeros(gridnz,gridtime,gridny);
compW = zeros(gridnz,gridtime,gridny);

for i = 1:gridtime
    SqueezeCompU = squeeze(windfield.u(:,i,:));
    compU(:,i,:) = flipud(SqueezeCompU');
    
    SqueezeCompV = squeeze(windfield.v(:,i,:));
    compV(:,i,:) = flipud(SqueezeCompV');
    
    SqueezeCompW = squeeze(windfield.w(:,i,:));
    compW(:,i,:) = flipud(SqueezeCompW');
end

clear windfield

%% Discretization

fullTime       = (dt*gridtime)-dt; %total time duration of the windfield
fullslicesTime = 0:dt:fullTime;
slicesDistance = fullslicesTime*Uref; % vector with distance between slices(m)
%Calculate the slice of each pattern point
icunt1 = 0;
for ipoint = 1:length(Y) % different time and slices for each point!!
    
    if input.flag_floating
        % finds the index of wind field slice where lidar measurements take place
        index = ((1+(icunt1*timeStep_Measurements)/dt):(size(compU,2)*(timestep_pat_vec/fullTime)):gridtime);
        % takes into account the position in x-direction due to floating dynamics
        slicesN     = floor(index + ((x_I{1,ipoint}(1,1:length(index))-ref_plane_dist)/Uref/dt));
        slicesN(slicesN<1) = NaN; % Slices outside the Grid in the beginning are NaN and discarded in next step
        slicesN(slicesN>gridtime) = NaN;% Slices outside the Grid in the End are NaN and discarded in next step
    else
        slicesN     = floor((1+(icunt1*timeStep_Measurements)/dt):(size(compU,2)*(timestep_pat_vec/fullTime)):gridtime); % Slices for each pattern point [index]. We use floor to avoid taking a slice that does not exist. Use high resolution WF!!!!
    end
    
    slicesTimeN = (icunt1*timeStep_Measurements:timestep_pat_vec:fullTime);% Slices for each pattern point [s]
    
    if length(slicesN) < length(floor((1+1:(size(compU,2)*(timestep_pat_vec/fullTime)):gridtime))) % check to make sure that the total length is correct compared to the total possible slices
        slicesN(end+1:length(floor((1+1:(size(compU,2)*(timestep_pat_vec/fullTime)):gridtime)))) = nan ;
    end
    if length(slicesTimeN)<length(0:timestep_pat_vec:fullTime) % check that the time vector has the same length for all points
        slicesTimeN(end+1:length(0:timestep_pat_vec:fullTime)) = nan ;
    end
    if ipoint>1 && length(slicesN)<length(slices(ipoint-1,:)) % check to make sure that the total length is correct compared to the other points
        slicesN(end+1:length(slices(ipoint-1,:))) = nan ;
    end
    
    slices(ipoint,:)     =  slicesN; %#ok<*AGROW> requested slices for scanning
    slicesTime(ipoint,:) =  slicesTimeN; % time of slices
    icunt1 = icunt1 +1;
    clear slicesN slicesTimeN
end

% Delete rows with NaNs in it. Follows the same logic as for non floating
if input.flag_floating
    [~ ,col] = find(isnan(slices));
    delete_slices = unique(col);
    slices(:,delete_slices) = [];
    slicesTime(:,delete_slices) = [];
    
    for i= 1:Pattern_n
        x_I{1,i}(:,delete_slices) = [];
        y_I{1,i}(:,delete_slices) = [];
        z_I{1,i}(:,delete_slices) = [];
        
        y_posL_I{1,i}(:,delete_slices) = [];
        z_posL_I{1,i}(:,delete_slices) = [];
        
        x_vel{1,i}(:,delete_slices) = [];
        y_vel{1,i}(:,delete_slices) = [];
        z_vel{1,i}(:,delete_slices) = [];
    end
else
    % Remove the patterns that include at least one nan. We keep only full pattern scans at the end
    SlicestoCut = [];
    NanMatSlices = (reshape(isnan(slices),size(slices,1),size(slices,2)));
    for ind_nan = 1:length(Y)  % get all the inices of nan values
        if any(NanMatSlices(ind_nan,:))
            SlicestoCut = [SlicestoCut,find(NanMatSlices(ind_nan,:) == 1)];
        end
    end
    % Cut slices and slicestime
    if ~isempty(SlicestoCut)
        MaxColumn  = size(slices,2)-min(SlicestoCut)+1;
        slices     = slices(:,1:end-MaxColumn);
        slicesTime = slicesTime(:,1:end-MaxColumn);
    end
end
distance_av_slice = distance_av_space/(distanceSlices); %transforming the distance in space to slices count

%% Rotate wind field for adding yaw and tilt inflow offsets

if input.yaw_angle ~= 0 && input.tilt_angle ~= 0
    [compU,compV,compW] = RotateWindfield(compU,compV,compW,input); %Function Not yet pushed
end

%% Extract points with time, spatial discretization and LOS

% floating trajectories
if input.flag_floating
    for num_tr = 1:length (Y)
        trajectory_fl(:,num_tr) = [{y_I{1,num_tr}};{z_I{1,num_tr}}]; %#ok<CCAT1,*SAGROW>  first line is y direction second is z. Origin is (Pos_LiDAR(1),Pos_LiDAR(2))==0
        trajectory_forAng_fl{1, num_tr} = trajectory_fl{1, num_tr}+y_posL_I{1,num_tr};
        trajectory_forAng_fl{2, num_tr} = trajectory_fl{2, num_tr}+z_posL_I{1,num_tr};
    end
end

% fixed trajectory is needed also in floating case for reconstruction
for num_tr = 1:length (Y)
    trajectory(:,num_tr) = [Y(num_tr);Z(num_tr)]; %#ok<*SAGROW>  first line is y direction second is z. Origin is (Pos_LiDAR(1),Pos_LiDAR(2))==0
    trajectory_forAng(:,num_tr) = [Y(num_tr)+Pos_LiDAR(1);Z(num_tr)+Pos_LiDAR(2)]; %#ok<*SAGROW>  Offset is included here to calculate the changed LOS. first line is y direction second is z. Origin is (Pos_LiDAR(1),Pos_LiDAR(2))==0
end

% loop over trajectory to find LOS angles (constant for all repeated measurement points)
if input.flag_apply_LOS == 1
    for i_ang = 1:length (Y)
        if input.flag_floating
            % floating angles
            angley_fl(i_ang,:) = atand(trajectory_forAng_fl{1,i_ang}/ref_plane_dist) ;
            anglez_fl(i_ang,:) = atand(trajectory_forAng_fl{2,i_ang}/ref_plane_dist) ;
        end
        % fixed angles
        angley(i_ang) = atand(trajectory_forAng(1,i_ang)/ref_plane_dist) ;
        anglez(i_ang) = atand(trajectory_forAng(2,i_ang)/ref_plane_dist) ;
    end
else
    for i_ang = 1:length (Y)
        if input.flag_floating
            angley_fl(i_ang,:) = trajectory_forAng_fl{1,i_ang};
            anglez_fl(i_ang,:) = trajectory_forAng_fl{2,i_ang};
            angley_fl(i_ang,:) = 0;
            anglez_fl(i_ang,:) = 0;
        end
        angley(i_ang) = 0;
        anglez(i_ang) = 0;
    end
end

%get the transformation matrices for each pattern point. look at https://en.wikipedia.org/wiki/Rotation_matrix
for iTra = 1:length(Y)
    %floating transformtion matrices
    if input.flag_floating
        for i_point=1:length(angley_fl)
            In_2_LOS_matrix_fl{iTra,i_point} = [cosd(anglez_fl(iTra,i_point))*cosd(angley_fl(iTra,i_point))   -sind(anglez_fl(iTra,i_point))   cosd(anglez_fl(iTra,i_point))*sind(angley_fl(iTra,i_point));
                                                sind(anglez_fl(iTra,i_point))*cosd(angley_fl(iTra,i_point))    cosd(anglez_fl(iTra,i_point))   sind(anglez_fl(iTra,i_point))*sind(angley_fl(iTra,i_point)) ;
                                                               -sind(angley_fl(iTra,i_point))                             0                     cosd(angley_fl(iTra,i_point))];
        end 
    end
    % %get the transformation matrices for each pattern point. look at https://en.wikipedia.org/wiki/Rotation_matrix
    In_2_LOS_matrix{iTra} = [cosd(anglez(iTra))*cosd(angley(iTra))   -sind(anglez(iTra))   cosd(anglez(iTra))*sind(angley(iTra));
                             sind(anglez(iTra))*cosd(angley(iTra))    cosd(anglez(iTra))   sind(anglez(iTra))*sind(angley(iTra)) ;
                                       -sind(angley(iTra))                     0                     cosd(angley(iTra))];
    LOS_2_In_matrix{iTra} =  In_2_LOS_matrix{iTra}^-1; % Inverse transformation: LOS_CS to Inertial_CS for reconstruction
end

if distance_av_slice ~= 0
    SliceVecInt = round((-distance_av_slice:distance_av_slice/points_av_slice:distance_av_slice))*distanceSlices;
    focus_distances = SliceVecInt+ref_plane_dist;
else
    focus_distances = ref_plane_dist;  % no slices to be averaged, single point measurement
end

%loop over planes to get points
for ifDist = 1:length(focus_distances)
    iplane = focus_distances(ifDist); % requested planes are all the planes along the distance (slicesDistance)
    
    for iTraj = 1:size(trajectory,2)
        if input.flag_apply_LOS ==1
            if input.flag_floating
                plane_traj{iTraj,ifDist}(1,:) = iplane*tand(angley_fl(iTraj,:))+ y_posL_I{1,iTraj}; % The position is added here as an offset to consider it in the coordinates of query points
                plane_traj{iTraj,ifDist}(2,:) = iplane*tand(anglez_fl(iTraj,:))+ z_posL_I{1,iTraj};
            else
                plane_traj{ifDist}(1,iTraj) = iplane*tand(angley(iTraj));
                plane_traj{ifDist}(2,iTraj) = iplane*tand(anglez(iTraj));  % this variable saves Y aand z points according to the plane
            end
        else   %dont move the trajectory projection if you dont have LOS
            if input.flag_floating
                plane_traj{iTraj,ifDist}(1,:) = y_I{1,iTraj};
                plane_traj{iTraj,ifDist}(2,:) = z_I{1,iTraj};
            else
                plane_traj{ifDist}(1,iTraj) = Y(iTraj);
                plane_traj{ifDist}(2,iTraj) = Z(iTraj); % this variable saves Y aand z points according to the plane
            end
        end
    end
end

LOS_points.slicesAv = round(-distance_av_slice:distance_av_slice/points_av_slice:distance_av_slice);

%Check for the case when we request only 1 slice to be averaged
if [isempty(LOS_points.slicesAv) || any(isnan(LOS_points.slicesAv))]  && distance_av_slice==0 %#ok<*NBRAK,BDSCA>
    LOS_points.slicesAv = 0;
end

for num_tr3 = 1:length (Y)
    LOS_points.slices(num_tr3,:)     = slices(num_tr3,:);
    LOS_points.slicesTime(num_tr3,:) = slicesTime(num_tr3,:);
    
    for iTraj2 = 1:size(plane_traj,2)
        for i_slice = 1:length(slices)
            if input.flag_floating
                LOS_points.Coor{num_tr3,i_slice}(:,iTraj2) = plane_traj{num_tr3,iTraj2}(:,i_slice); %this variable saves coordinates according to trajectory points
            else
                LOS_points.Coor{num_tr3,i_slice}(:,iTraj2) = plane_traj{1,iTraj2}(:,num_tr3);
            end
        end
    end
end

% extract the measured slices as well as the full time series. NO LOS here
% The output has a vector of values for each measured point representing 
% the inertial values of the wind speed in the inertial frame

[VFinalTotal_U,VFinalTotal_Time_U,~,~] = interpolationFun(compU,LOS_points,gridy,gridz,fullTime,dt,type_interpolation_2, Y, Z);
[VFinalTotal_V,VFinalTotal_Time_V,~,~] = interpolationFun(compV,LOS_points,gridy,gridz,fullTime,dt,type_interpolation_2, Y, Z);
[VFinalTotal_W,VFinalTotal_Time_W,~,~] = interpolationFun(compW,LOS_points,gridy,gridz,fullTime,dt,type_interpolation_2, Y, Z);

%add translational velocity for floating case
if input.flag_floating
    for num_tr = 1:length (Y)
        VFinalTotal_Time_U{1,num_tr} = VFinalTotal_Time_U{1,num_tr}-repmat(x_vel{1,num_tr},length(LOS_points.slicesAv),1);
        VFinalTotal_Time_V{1,num_tr} = VFinalTotal_Time_V{1,num_tr}-repmat(y_vel{1,num_tr},length(LOS_points.slicesAv),1);
        VFinalTotal_Time_W{1,num_tr} = VFinalTotal_Time_W{1,num_tr}-repmat(z_vel{1,num_tr},length(LOS_points.slicesAv),1);
    end
end

%% LOS transformations for all cases

for ind_LOS = 1:length(Y) % loop over the beam
    for ind_slice = 1:length(VFinalTotal_Time_U{ind_LOS}) % loop over the measurements in time e.g 200  
        for P_av= 1:length(LOS_points.slicesAv) % loop over all points along probe volume LOS e.g. 20
            if input.flag_floating
                %multiplyng all the slice with the transformation matrix
                VFinalTotal_Time_LOS_vec{ind_LOS, ind_slice}(:,P_av) =  In_2_LOS_matrix_fl{ind_LOS,ind_slice} * ...
                    [VFinalTotal_Time_U{ind_LOS}(P_av,ind_slice);VFinalTotal_Time_V{ind_LOS}(P_av,ind_slice);VFinalTotal_Time_W{ind_LOS}(P_av,ind_slice)];
                                                                                                     
            else
                %multiplyng all the slice with the transformation matrix
                VFinalTotal_Time_LOS_vec{ind_LOS, ind_slice}(:,P_av) =  In_2_LOS_matrix{ind_LOS} * ...
                    [VFinalTotal_Time_U{ind_LOS}(P_av,ind_slice);VFinalTotal_Time_V{ind_LOS}(P_av,ind_slice);VFinalTotal_Time_W{ind_LOS}(P_av,ind_slice)];
            end                                                                                                      
        end
    end
end
% VFinalTotal_Time_U includes: Each cell is row1:u row2:v, row2:w coulumns
% are points in the probe volume. The rows of the struct are the beams and
% colums are the measurements in time. No probe volume averaging until here


%% Apply probe volume weighting                                                                                                                                                                            
[VFinalTotal_Time_LOS_vec] = Weighting_Probe(VFinalTotal_Time_LOS_vec, Y); % here VFinalTotal_Time_U changes so that each cell has only 1 column 
%using the average. If the raw scatter is required has to be retrieved here.  

%% Application of noise to the measured points. This should be a function
% noise
if input.flag_apply_noise == 1
    for ind_LOS= 1:length(Y)
           LOS_vec= cell2mat(VFinalTotal_Time_LOS_vec(ind_LOS,:)); % create each beam individually in one variable
           LOS_vec(1,:) = awgn(LOS_vec(1,:),noise_U);
           LOS_vec(2,:) = awgn(LOS_vec(2,:),noise_V);
           LOS_vec(3,:) = awgn(LOS_vec(3,:),noise_W);
        for ind_slice=1:size(VFinalTotal_Time_LOS_vec,2) % converting back to the cell format as before. Could be done faster..
            VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(1,1) = LOS_vec(1,ind_slice);
            VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(2,1) = LOS_vec(2,ind_slice);
            VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(3,1) = LOS_vec(3,ind_slice);
        end
    end
end

%% decompose to components. Not use at the moment but useful for future applicaitons
for ind_LOS = 1:length(Y)
    for ind_slice = 1:length(VFinalTotal_Time_U{ind_LOS})
        VFinalTotal_Time_LOS_U{ind_LOS}(ind_slice) = VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(1); %#ok<NASGU>
        VFinalTotal_Time_LOS_V{ind_LOS}(ind_slice) = VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(2); %#ok<NASGU>
        VFinalTotal_Time_LOS_W{ind_LOS}(ind_slice) = VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(3); %#ok<NASGU>
    end
end

%% Wind Field Reconstruction
[VFinalTotal_Time_U, VFinalTotal_Time_V, VFinalTotal_Time_W] = WindFieldReconstruction(VFinalTotal_Time_LOS_vec, LOS_2_In_matrix, input, Y);


%% Calculate REWS

% distances of all points in the turbsim grid from center assuming hub center =0
for I = 1:gridny
    for II = 1:gridnz
        Distances_of_Points_In_Plane(II,I) = sqrt((gridy(I)).^2+(gridz(II)).^2); %Matrix of distances from center of rotor to each point in the grid
    end
end

% distances of all points of the lidar pattern from center assuming hub center =0
for ind_length = 1:length(Y)
    Distances_LiDAR_Points(ind_length) = sqrt(Y(ind_length).^2+(Z(ind_length)).^2); %Matrix of distances from center of rotor to each point in the grid
end

if input.flag_apply_weightREWS == 1
    dist_REWS = rotor_radius*[dist_REWS_nd 10]; %convert ND spanwise to meters and treat the point outside with 0 weight
    Wi        = [Wi 0];
    % calculate weigths for all the turbsim  grid points
    for I = 1:gridny
        for II = 1:gridnz
            Weights_of_Points_In_Plane(II,I) = interp1(dist_REWS,Wi,Distances_of_Points_In_Plane(II,I)); %Matrix of weights for all grid points
        end
    end
    weigthTot_grid = sum(sum(Weights_of_Points_In_Plane)); %sum of weights needed for weighted average
    
    % calculate weigths for all lidar pattern points
    for ind_length = 1:length(Y)
        Weights_of_Lidar_In_Plane(ind_length) = interp1(dist_REWS,Wi,Distances_LiDAR_Points(ind_length));%Matrix of weights for all lidar pattern points
    end
    weigthTot_lidar = sum(sum(Weights_of_Lidar_In_Plane)); %sum of weights needed for weighted average
end

%%%%%%%%%%%%%%%%%%%%%%%%% For FullWindField data: %%%%%%%%%%%%%%%%%%%%%%%%%

%Transform Slices
rad_values = Distances_of_Points_In_Plane<=rotor_radius;       %keep only points inside the rotor
for ind_slicesDistance = 1:length (slicesDistance) % loop over all the slices
    % Take complete slice from turbsim:
    ExSlice_U = squeeze(compU(:,ind_slicesDistance,:)); % Values of the selected slice CompU
    ExSlice_U = ExSlice_U.*rad_values; %remove points out of rotor
    if input.flag_apply_weightREWS == 1
        ExSlice_U = ExSlice_U.*Weights_of_Points_In_Plane;         %multiply with weights
        REWS.fullWF.TS(ind_slicesDistance) = sum(sum(ExSlice_U))/weigthTot_grid;
    else
        noZeroSlice = nonzeros(ExSlice_U);
        REWS.fullWF.TS(ind_slicesDistance) = mean(noZeroSlice,'omitnan');
    end
end
REWS.fullWF.mean = mean(REWS.fullWF.TS,'omitnan');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% For the pattern of LiDAR points %%%%%%%%%%%%%%%%%%%%%%%%%
rad_valuesLiDAR = Distances_LiDAR_Points <= rotor_radius; %logical index of pattern points inside the rotor
for ind_LiDAR_REWS = 1:size(slicesTime,2) % here there is something to fix:  pointsToAverage sometimes (when TstepPattern=time step of the Original WF) gives errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ind_point_LiDAR_REWS = 1:length(Y)
        pointsToAverage(:,ind_point_LiDAR_REWS) = VFinalTotal_Time_U{ind_point_LiDAR_REWS}(1,ind_LiDAR_REWS); % we take points of each Time series in the correspondent slice and mean it.
    end
    pointsToAverage = rad_valuesLiDAR.*pointsToAverage; % remove points out of the rotor
    if input.flag_apply_weightREWS == 1
        pointsToAverage = Weights_of_Lidar_In_Plane.*pointsToAverage;
        REWS.lidar.TS(ind_LiDAR_REWS) = sum(pointsToAverage)/weigthTot_lidar;
    else
        REWS.lidar.TS(ind_LiDAR_REWS) = mean(nonzeros(pointsToAverage),'omitnan'); % removing points that are outside the rotor area
    end
end
REWS.lidar.mean = mean(REWS.lidar.TS,'omitnan');
if timeStep_Measurements ~= 0
    REWS.lidar.TSTime = min(min(slicesTime))+ timestep_pat_vec/2:timestep_pat_vec:max(max(slicesTime))+timeStep_Measurements; % We suppose that measurements are in the center of the pattern (half the time of pattern time step). If do not create the same number of points in time we add the points
    REWS.lidar.TSTime = interp1(fullslicesTime,fullslicesTime,REWS.lidar.TSTime,'nearest');
    if length(REWS.lidar.TSTime)<length(slicesTime) && length(REWS.lidar.TSTime)<length(REWS.lidar.TS) % Consider to clean this up and merge with the floating case above
        No_difference_points = length(slicesTime)-length(REWS.lidar.TSTime);
        for idifpoint = 1:No_difference_points
            REWS.lidar.TSTime = [REWS.lidar.TSTime REWS.lidar.TSTime(end)+ timestep_pat_vec]; % 
        end
    end
else   % when we have all beams synchronized we dont need to average time
    REWS.lidar.TSTime = slicesTime(1,:);
end

%% Calculate Shear power law exponent from measurements and full windfield
%TODO: maybe add an option to calculate shear on every nth slice of the full field to reduce time

zero_valueZ = find(gridz==0);
if isempty(zero_valueZ)
    [~,zero_valueZ] = min(abs(gridz));   % dirty fix for case where the grid does not include (0,0)
end

z_vec_Shear = (gridz)+Zh; %create vector of heights
Vhub_shear = zeros(gridtime,1); % pre assign variable
% calculate power law for full turbsim data:
for ind_sliceLaW = 1:gridtime  %for slices
    
    V_slice = squeeze(compU(:,ind_sliceLaW,:)); % velocity of point in slices
    Vhub_shear(ind_sliceLaW)  = mean(V_slice(zero_valueZ,:));
    v_hor   = mean(V_slice,2);     % take the average of all points in each horizontal line
    
    % find least square fit for the average vertical line
    fcn = @(alphaPL) sum((Vhub_shear(ind_sliceLaW)*(z_vec_Shear/Zh).^(alphaPL) - v_hor').^2); % least square defintion f
    [s,~,~] = fminsearch(fcn,0.14);        % Minimise Least-Squares error
    ShearPL.fullWF.TS(ind_sliceLaW) = s; %#ok<*SAGROW>
    
end
ShearPL.fullWF.Mean = mean(ShearPL.fullWF.TS); % total shear of the wind field

%%%%%%% Calculate power law exponent from LiDAR measurements.%%%%%%%%%%%%%%

% We suppose that each pattern can be projected to one plane in order to calculate shear. This plane is the middle of the scan pattern in time
Match_index = FindSameValuesAndAverageV(Z);   % Calculate index of repeated heights in order to calculate mean of each height for calculating EXPONENT LAW
%     Z_lid_point = ceil(length(LOS_points.slicesAv)/2); %find the middle distance of the multiple ranges  where we project all the velocities
[~,Z_sorted_ind] = sort(Z);
Match_index_mat  = cell2mat(Match_index);
Match_index2 = Match_index;
for indSlice = 1:size(VFinalTotal_Time_U{1},2)
    for ind_LOS = 1:length(Z)
        iV_vec_lidar(ind_LOS) =  VFinalTotal_Time_U{ind_LOS}(indSlice);
    end
    
    % create vector of heights and speeds sorted by speeds from lowerto higher
    % TODO: THIS PART IS SPAGHETTI CODE WE HAVE TO RECONSIDER FOR A BETTER SOLUTION
    for iZ = 1:length(Z)
        iZint = Z_sorted_ind(iZ);                      % take the index of the original height vector corresponding to the order of sorted heights
        if any(Match_index_mat == iZint)                 % check if the height is included in the matching indices
            for ind_LAW_NO_LOS=1:size(Match_index,2)   % loop ove all matching groups of heights
                iMatchIndex = Match_index2{ind_LAW_NO_LOS};
                if any(iMatchIndex == iZint) % check if my value is in this group of matching values
                    for ind_LAW_NONLOS2 = 1:length(iMatchIndex) %loop over the matching elements
                        if ind_LAW_NONLOS2 == 1 % if it is the first of the matching group assign the mean value and the height in the second sorted index
                            Z1_sh(iZ) = Z(iZint);
                            V1_sh(iZ) = mean(iV_vec_lidar(iMatchIndex),'omitnan');
                        else
                            Match_index2{ind_LAW_NO_LOS}(ind_LAW_NONLOS2) = nan; %if it not the first just clear the value to avoid repeateing it next
                        end
                    end
                end
            end
            if length(Z1_sh) < iZ
                Z1_sh(iZ) = nan; %if it passed all htes tests and it is nothere assign nan to keep length
                V1_sh(iZ) = nan;
            end
            
        else % if it is not in the matching groups assign to a new variable the values height and velocity according to the sorted values
            Z1_sh(iZ) = Z(iZint);
            V1_sh(iZ) = iV_vec_lidar (iZint);
        end
    end
    Vshear_final = V1_sh(~isnan(V1_sh));
    Zshear_final = Zh + Z1_sh(~isnan(Z1_sh));
    z_vertical_min=min(abs(Zshear_final-Zh)); % finds the minimum vertical distance from hub height
    z_vertical=abs(Zshear_final-Zh); % vector of distances from hub
    
    [~,idx_ZShear_min] = find(z_vertical==z_vertical_min); % finds index of minimum vertical distance from hub height
    Vhub_shear = mean(Vshear_final(idx_ZShear_min)); % finds the mean of all points with minimum distance
    clear V1_sh Z1_sh
    
    % Calculate the approximate Vhub (or take it from turbsim...)
    % Find least square fit for the average vertical line
    if length(Vshear_final) == length(Zshear_final)
        fcn     = @(alphaPL2) sum(( Vhub_shear.*(Zshear_final/Zh).^(alphaPL2) - Vshear_final).^2); % least square defintion f
        [s,~,~] = fminsearch(fcn, 0.14);        % Minimise Least-Squares error
        ShearPL.lidar.TS(:,indSlice) = s;       %#ok<*SAGROW>
    else
        ShearPL.lidar.TS(:,indSlice) = nan;     %#ok<*SAGROW> % in case there are many nans there is not enough data for a height so discard measurement
    end
    Match_index2 = Match_index;
end
ShearPL.lidar.Mean = mean(ShearPL.lidar.TS,'omitnan'); % total shear of the wind field
if timeStep_Measurements ~= 0
    ShearPL.lidar.TSTime = min(min(slicesTime))+ timestep_pat_vec/2:timestep_pat_vec:max(max(slicesTime))+timeStep_Measurements; % We suppose that measurements are in the center of the pattern (half the time of pattern time step). If do not create the same number of points in time we add the points
    ShearPL.lidar.TSTime = interp1(fullslicesTime,fullslicesTime,ShearPL.lidar.TSTime,'nearest');
    if length(ShearPL.lidar.TSTime)<length(slicesTime) && length(ShearPL.lidar.TSTime)<length(ShearPL.lidar.TS)
        No_difference_points = length(slicesTime)-length(ShearPL.lidar.TSTime);
        for idifpoint = 1:No_difference_points
            ShearPL.lidar.TSTime = [ShearPL.lidar.TSTime ShearPL.lidar.TSTime(end)+ timestep_pat_vec]; % 
        end
    end
else   % when we have all beams synchronized we dont need to average time
    ShearPL.lidar.TSTime = slicesTime(1,:);
end

%% Statistics
% obtain metrics of convergence between original and constrained WindField, mean of
% error, STDV of the current Time series and STDV of the
% error...

[statisticsOut.U] = statisticsFun(Y,VFinalTotal_U,VFinalTotal_Time_U,slices,slicesTime);
[statisticsOut.V] = statisticsFun(Y,VFinalTotal_V,VFinalTotal_Time_V,slices,slicesTime);
[statisticsOut.W] = statisticsFun(Y,VFinalTotal_W,VFinalTotal_Time_W,slices,slicesTime);

%% Resample Data (check the missmatch in time at the end of the series...!!)

if input.flag_resampling == 1
    slicesTime2  = slicesTime;
    clear slicesTime
    for iPat = 1:length(Y)
        %CompU
        VFinalTotal_Time_U{iPat} = interpft(VFinalTotal_Time_U{iPat},length(VFinalTotal_Time_U{iPat})*resampling_factor);
        VFinalTotal_Time_U{iPat} = VFinalTotal_Time_U{iPat}(1:end-resampling_factor+1);
        %CompV
        VFinalTotal_Time_V{iPat} = interpft(VFinalTotal_Time_V{iPat},length(VFinalTotal_Time_V{iPat})*resampling_factor);
        VFinalTotal_Time_V{iPat} = VFinalTotal_Time_V{iPat}(1:end-resampling_factor+1);
        %CompW
        VFinalTotal_Time_W{iPat} = interpft(VFinalTotal_Time_W{iPat},length(VFinalTotal_Time_W{iPat})*resampling_factor);
        VFinalTotal_Time_W{iPat} = VFinalTotal_Time_W{iPat}(1:end-resampling_factor+1);
        
        slicesTime(iPat,:)       = slicesTime2(iPat,1):(slicesTime2(iPat,2)-slicesTime2(iPat,1))/resampling_factor:slicesTime2(iPat,end);
    end
    % Matching the time vector with TurbSim time series:
    slicesTime = slicesTime(:,1:length(VFinalTotal_Time_U{1}));
    timestep_pat_vec = timestep_pat_vec/resampling_factor;
else
    [~,colNaNslicesTime] = find(isnan(slicesTime));
    ncolsdelete = unique(colNaNslicesTime);
    slicesTime  = slicesTime(:,1:(end-length(ncolsdelete)));
end

%% Create and save .mat lidar output

Output.REWS  = REWS;
Output.Shear = ShearPL;

if input.flag_apply_noise == 1
    Output.Parameter.Noise  = [noise_U ;noise_V; noise_W];
else
    Output.Parameter.Noise  = [] ;
end
Output.statistics       = statisticsOut;
Output.TS.fullWF.time   = fullslicesTime;
Output.TS.fullWF.duration = fullslicesTime(end)+dt; % This is the length of the full wind field in seconds. dt is added to start the timeseries from zero
Output.Pattern.Coord    = [Y;Z];
Output.Pattern.refplane = ref_plane_dist; %like focus distance
Output.Pattern.timestep_pat_vec  = timestep_pat_vec;
Output.Pattern.timeStep_Measurements = timeStep_Measurements;
Output.Pattern.distance_av_slice = distance_av_slice;
Output.Pattern.points_av_slice   = points_av_slice;
Output.Pattern.timestep_pat_vec  = timestep_pat_vec;
Output.Pattern.name              = input.PatternNames{curFileInfo.values{find(strcmp(curFileInfo.variables{1, 1},'Pat'))}};
Output.Parameter.rotor_radius = rotor_radius; % Radius of the Rotor [m]
Output.Parameter.Pos_Lidar    = input.Pos_LiDAR;

for iPat= 1:length(Y)
    Output.TS.fullWF.Uval{iPat} = VFinalTotal_U{iPat};
    Output.TS.fullWF.Vval{iPat} = VFinalTotal_V{iPat};
    Output.TS.fullWF.Wval{iPat} = VFinalTotal_W{iPat};
    Output.TS.lidar.Uval{iPat}  = VFinalTotal_Time_U{iPat};
    Output.TS.lidar.Vval{iPat}  = VFinalTotal_Time_V{iPat};
    Output.TS.lidar.Wval{iPat}  = VFinalTotal_Time_W{iPat};
    Output.TS.lidar.time{iPat}  = slicesTime(iPat,:);
end
Output.TS.fullWF.nGridY = gridny;
Output.TS.fullWF.nGridZ = gridnz;
Output.TS.fullWF.dy = dy;
Output.TS.fullWF.dz = dz;

%save
save_data_full_path = [input.LidarOutput_dir curFileInfo.name '.mat'];
save(save_data_full_path,'Output')
