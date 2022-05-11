%% Header
%
% Create the .TimeSer input file for turbsim constrained turbulence
%
% TODO:
% 1) find a way to cut the time series to avoid the Delta frequency error
% from turbsim
% 2) Consider move the extension of Vinput earlier in the code to make sure
% all time series start from zero, using NaN.
%
%
% V.Pettas/F.Costa/M.Gräfe
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function [VFinal,t_lidar,Y1,Z1] = Create_Turbsim_Vinput(LidarOutput,input)
%First sort with heigths
Z = LidarOutput.Pattern.Coord (2,:);
Y = LidarOutput.Pattern.Coord (1,:);

U_velTS_init = LidarOutput.TS.lidar.Uval ;
V_velTS_init = LidarOutput.TS.lidar.Vval ;
W_velTS_init = LidarOutput.TS.lidar.Wval ;
t_lidar = LidarOutput.TS.lidar.time{1, 1};
nComp   = input.nComp ;

% DO SOMETHING FOR THE CASE THAT POINTS ARE NOT THE SAME time stamp
[Z1,SortedIndex] = sort(Z);
Y1 = Y(SortedIndex);
for ind_sort = 1:length(Y1)
    U_velTS{ind_sort} = U_velTS_init{SortedIndex(ind_sort)};
    V_velTS{ind_sort} = V_velTS_init{SortedIndex(ind_sort)}; %#ok<*AGROW>
    W_velTS{ind_sort} = W_velTS_init{SortedIndex(ind_sort)};
end
Z1 = Z1+input.Zh ; % adding hub height

% Take the Time Series depending on the components requested
switch nComp
    case 1
        VFinal = [];
        for iTurb = 1:length(Y)
            VFinal(:,iTurb) = (U_velTS{iTurb}');
        end
        nComp_TurbSim = nComp; %#ok<*NASGU>
    case 2
        for iTurb = 1:length(Y)
            VFinal2{iTurb} = horzcat(U_velTS{iTurb}',(V_velTS{iTurb})');
        end
        VFinal = [];
        for iTurb2 = 1:length(Y)
            VFinal = [VFinal,(VFinal2{iTurb2})];
        end
        nComp_TurbSim = nComp;
    case 3
        for iTurb = 1:length(Y)
            VFinal3{iTurb} = horzcat(U_velTS{iTurb}',(V_velTS{iTurb})',(W_velTS{iTurb})');
        end
        VFinal = [];
        for iTurb2 = 1:length(Y)
            VFinal = [VFinal,(VFinal3{iTurb2})];
        end
        nComp_TurbSim = nComp;
end

%     check for Nans and remove them if turbsim cannot handle them
for iRowTS = 1:size(VFinal,1)
    if any(isnan(VFinal(iRowTS,:)))
        CntNan(iRowTS) = 0;
    else
        CntNan(iRowTS)=1;
    end
end

VFinal = (VFinal(find(CntNan==1),:));  %#ok<FNDSB>

% check if first timestep is 0, otherwise extend V input because Turbsim requires to have the first entry to be at 0s
firstTS = t_lidar(1,1);
delta = t_lidar(1,2)-t_lidar(1,1); % get sampling rate
if t_lidar(1,1)>0 % fill lines until first entry is zero, first existing value is repeated
    fillline = ones(length(0:delta:firstTS-delta),length(Y)*nComp).*VFinal(1,:);
    VFinal = [ fillline; VFinal]; % extent values
    t_lidar = [ 0:delta:firstTS-delta, t_lidar]; % extent timestamps
end
clearvars fillline

% extend to last valid timestamp required by Turbsim by adding the last
% available value for each beam
last_timestep = LidarOutput.TS.fullWF.duration-delta;
if t_lidar(1,end)<last_timestep % check if the last timestep of time series is the last entry required by turbsim (AnalysisTime-delta)
    fillline = ones(length( t_lidar(1,end)+delta:delta:last_timestep),length(Y)*nComp).*VFinal(end,:); % fill missing lines with mean value
    VFinal = [VFinal; fillline ]; % extent values
    t_lidar = [t_lidar, t_lidar(1,end)+delta:delta:last_timestep]; % extent timestamps
end

% TODO:
%     % Call function for generating turbsim timeusr inputs
%     % cut the time series to matcvh with the original:
%     %Here have to add the lines of code that finds the time in which all
%     %time series coincide: (script PruebaTotalTimeSeries in the desktop)
%
%     %     time_row=find(t==Analysistime);
%         %     Analysistime=298; % I added this to match the time series vector. Have to figure out how to do it automatically. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     writeTurbSimTimeSeriesInput(input_file, windfieldfile, VFinal(1:time_row,:), t_lidar(1:time_row,:), Y1, Z1, RefNode,nComp_TurbSim); % Write time series for TurbSim
%     %     disp(strcat('Time series of the wind field  ', windfieldfile(1:end-4),' for TurbSim is done'))
%     %Take the constraining timeseries, pass it to Turbsim and create the
%     %template for turbsim
%     GridHt=gridnz*dz-dz;
%     GridWt=gridny*dy-dz;
%
%     % writing .inp file for turbsim
%     processingTurbsim(TSnam2,SaveDirectory_TurbSim,OutSaveName,RandSeed,gridnz,gridny,timestep_pat_vec,Analysistime,UsableTime,Zh,GridHt,GridWt,IECstandard,statisticsOut.U.lidar.TI_mean,RefHt,statisticsOut.U.lidar.Mean_TS_LiDAR,ShearPL.lidar.Mean,WindProfileType,SCMod1)
%
%     disp([windfieldfile(1:end-4) ' template of the time series is done.' ])
