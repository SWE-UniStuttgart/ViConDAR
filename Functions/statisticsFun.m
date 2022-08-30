%% Header
%
% Calculate all statistics from the lidar measured time series and original wind
% field to store them in the Output variable. This function can be modified to
% introduce new desired outputs.
%
% V.Pettas/F.Costa/M.Graefe
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019


function [statisticsFunOut] = statisticsFun(Y,VFinalTotal_full,VFinalTotal_Lidar,slices,slicesTime)

for indstat=1:length(Y)
    %Mean
    MeanTS(indstat,:)       = mean(VFinalTotal_full{indstat},'omitnan');      %mean Of each full Time Series
    MeanTS_LiDAR(indstat,:) = mean(VFinalTotal_Lidar{indstat},'omitnan'); %mean Of each LiDAR Time Series
    
    %Error
    Error{indstat}     = (VFinalTotal_full{indstat} (slices(indstat,:))- VFinalTotal_Lidar{indstat}); %#ok<*AGROW> %Error between TurbSim-Ts and each point in the Pattern-TS:
    RMSE(indstat)      = sqrt(sum( (VFinalTotal_full{indstat}(slices(indstat,:)) - VFinalTotal_Lidar{indstat}).^2) / length(VFinalTotal_Lidar{indstat}) ); % RMSE for each TS of the pattern
    LidarTime{indstat} = slicesTime(indstat,:);
    
    ABSOLUT_Error{indstat}  = abs(Error{indstat}); %Absolut error 
    MeanError(indstat)      = mean(Error{indstat}); % mean error
    MeanAbsError(indstat)   = mean(ABSOLUT_Error{indstat});
    Max_Error(indstat)      = max(Error{indstat}); % Maximum of the error
    Min_Error(indstat)      = min(Error{indstat});% Minimum of the error
    Variance_Error(indstat) = var(Error{indstat});
    STDV_Error(indstat)     = std(Error{indstat});
    STDV_LiDAR(indstat)     = std(VFinalTotal_Lidar{indstat},'omitnan'); %STDV LiDAR Pattern-TS
    STDV_TS(indstat)        = std(VFinalTotal_full{indstat},'omitnan'); % STDV TurbSim-TS
    VAR_LiDAR(indstat)      = var(VFinalTotal_Lidar{indstat},'omitnan'); %VARIANCE LiDAR Pattern-TS
    VAR_TS(indstat)         = var(VFinalTotal_full{indstat},'omitnan'); % VARAIANCE TurbSim-TS
    TI_LiDAR(indstat)       = STDV_LiDAR(indstat) / MeanTS_LiDAR(indstat,:)    ;
    TI_TS(indstat)          = STDV_TS(indstat) / MeanTS(indstat,:)    ; %#ok<*NASGU>
end
MeanTotal_Full_TS  = mean(MeanTS); %Mean of the Full TimeSeries pattern
MeanTotal_LiDAR_TS = mean(MeanTS_LiDAR); %Mean of the LiDAR pattern
Mean_STDV_LiDAR    = mean(STDV_LiDAR);% Mean STDV LiDAR pattern
Mean_STDV_Full_TS  = mean(STDV_TS);% Mean STDV Full TimeSeries pattern
TI_mean_LiDAR_TS   = Mean_STDV_LiDAR/MeanTotal_LiDAR_TS; % TI LiDAR pattern
TI_mean_Full_TS    = Mean_STDV_Full_TS/MeanTotal_Full_TS;  % TI Full TimeSeries

%% Outputs

% Error
statisticsFunOut.error.Error          = MeanError;
statisticsFunOut.error.Error_TS       = Error;
statisticsFunOut.error.Abs_Error      = ABSOLUT_Error;
statisticsFunOut.error.Mean_Abs_Error = MeanAbsError;
statisticsFunOut.error.Max_error      = Max_Error;
statisticsFunOut.error.Min_Error      = Min_Error;
statisticsFunOut.error.STDV_Error     = STDV_Error;
statisticsFunOut.error.Variance_Error = Variance_Error;
statisticsFunOut.error.error_time     = LidarTime;

% LiDAR
statisticsFunOut.lidar.Mean_all_LiDAR    = MeanTotal_LiDAR_TS;
statisticsFunOut.lidar.Mean_Points_LiDAR = MeanTS_LiDAR';
statisticsFunOut.lidar.STDV_LiDAR        = STDV_LiDAR;
statisticsFunOut.lidar.VAR_LiDAR         = VAR_LiDAR;
statisticsFunOut.lidar.RMSE              = RMSE;
statisticsFunOut.lidar.TI                = TI_LiDAR;
statisticsFunOut.lidar.TI_mean           = TI_mean_LiDAR_TS;
statisticsFunOut.lidar.Mean_STDV         = Mean_STDV_LiDAR;

% FullWF
statisticsFunOut.fullWF.Mean_points_full = MeanTS';
statisticsFunOut.fullWF.Mean_all_full    = MeanTotal_Full_TS;
statisticsFunOut.fullWF.STDV_TS = STDV_TS;
statisticsFunOut.fullWF.VAR_TS  = VAR_TS;
statisticsFunOut.fullWF.TI      = TI_TS;
statisticsFunOut.fullWF.TI_mean = TI_mean_Full_TS;
statisticsFunOut.fullWF.Mean_STDV = Mean_STDV_Full_TS;

end