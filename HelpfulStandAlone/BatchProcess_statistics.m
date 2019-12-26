%% Header
%
% Dirty script to batch process all lidar measurements and create contours from
% statistics included in the lidar output files. It was done for a specific 
% application but can be generalized. 
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

%%
clear all %#ok<*CLALL>
clc
close all

%% IO definition
%Here we have to define the name pattern and all the expected variations

Base    = 'DTU10MW';
Shear   = 'Sh';
seed    = 'SD';
Vmean   = 'V';
TurbInt = 'TI';
PatternName = '9P_Cross';
TimestepPattern = 'Tp';
TimestepMeas    = 'Tm';

% lidar direcotrz
LidarDir = '../LidarOutput/';

%Choose which value you want to plot (index)
PlotNam = {'REWS mean of the TS abs error [%]' 'Shear mean of the TS abs error [%]' 'TI Error from mean of all points[%]' 'meanAbsError Error of the mean TS [%]' 'Error of the Mean Shear' 'Error of the mean REWS' 'Shear median of the TS abs error [%]' 'Shear mode of the TS abs error [%]'}; %corresponds to the outputs chosen see line 91-95

% We can combine in total 2 variations
Var1 = 5:5:25; % shear variation1
Var1Name = 'Shear';

Var2 = 5:5:30; % TI variation 2
Var3 = Var2/5;        % Tm variees only for names add the number of point in pattern here
Var2Name = 'TI';
varPat = 0; % flag to indicate that timestep of pattern is a variable (naming conventions)

% if you have multiple seeds you want to average:
SeedVar = 1:6;

%% Create the list of names

% if pattern timestep is a varuable usethis
if varPat== 1
    cnt = 1;
    for iVar1 = 1:length(Var1)
        for iVar2 = 1:length(Var2)
            for iSeed =SeedVar;
                if  mod(Var2(iVar2),1)==0 && mod(Var3(iVar2),1)==0
                    curNamInt = strrep([ Base '_' Shear num2str(Var1(iVar1),'%02.f') '_' seed num2str(SeedVar(iSeed),'%02.f') '_' Vmean '08_' TurbInt  '05_' PatternName '_' TimestepPattern num2str(Var2(iVar2),'%02.f') '_' TimestepMeas num2str(Var3(iVar2),'%02.f')],'.','d');
                elseif mod(Var2(iVar2),1)~=0 && mod(Var3(iVar2),1)==0
                    curNamInt = strrep([ Base '_' Shear num2str(Var1(iVar1),'%02.f') '_' seed num2str(SeedVar(iSeed),'%02.f') '_' Vmean '08_' TurbInt  '05_' PatternName '_' TimestepPattern num2str(Var2(iVar2),'%02.1f') '_' TimestepMeas num2str(Var3(iVar2),'%02.f')],'.','d');
                elseif mod(Var2(iVar2),1)==0 && mod(Var3(iVar2),1)~=0
                    curNamInt = strrep([ Base '_' Shear num2str(Var1(iVar1),'%02.f') '_' seed num2str(SeedVar(iSeed),'%02.f') '_' Vmean '08_' TurbInt  '05_' PatternName '_' TimestepPattern num2str(Var2(iVar2),'%02.f') '_' TimestepMeas num2str(Var3(iVar2),'%02.1f')],'.','d');
                elseif mod(Var2(iVar2),1)~=0 && mod(Var3(iVar2),1)~=0
                    curNamInt = strrep([ Base '_' Shear num2str(Var1(iVar1),'%02.f') '_' seed num2str(SeedVar(iSeed),'%02.f') '_' Vmean '08_' TurbInt  '05_' PatternName '_' TimestepPattern num2str(Var2(iVar2),'%02.1f') '_' TimestepMeas num2str(Var3(iVar2),'%02.1f')],'.','d');
                end
                ListLidar{cnt,1} = [LidarDir curNamInt '.mat'];
                ListLidar{cnt,2} = Var1(iVar1);
                ListLidar{cnt,3} = Var2(iVar2);
                cnt=cnt+1;
            end
        end
    end
    
else % if you dont have timestep of pattern as a variable use this
    cnt=1;
    for iVar1 = 1:length(Var1)
        for iVar2 = 1:length(Var2)
            for iSeed = SeedVar;
                curNamInt = strrep([ Base '_' Shear num2str(Var1(iVar1),'%02.f') '_' seed num2str(SeedVar(iSeed),'%02.f') '_' Vmean '08_' TurbInt  num2str(Var2(iVar2),'%02.f') '_' PatternName '_' TimestepPattern '2d3_' TimestepMeas '0d3'],'.','d'); 
                ListLidar{cnt,1} = [LidarDir curNamInt '.mat'];
                ListLidar{cnt,2} = Var1(iVar1);
                ListLidar{cnt,3} = Var2(iVar2);
                cnt=cnt+1;
            end
        end
    end
end

%% load and save all usefull values in one cell
for iNam = 1:size(ListLidar,1)
        try
    load(ListLidar{iNam,1}); % Output variable created REWS Shear TI meanAbsError
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ListLidar{iNam,4} = 100*mean(abs(Output.REWS.lidar.TS-Output.REWS.fullWF.TS(find(ismember(Output.TS.fullWF.time,Output.REWS.lidar.TSTime))))/Output.REWS.fullWF.mean); %#ok<*FNDSB,*SAGROW>  %Mean of the time series error of REWS
    ListLidar{iNam,5} = 100*mean(abs(Output.Shear.lidar.TS-Output.Shear.fullWF.TS(find(ismember(Output.TS.fullWF.time,Output.Shear.lidar.TSTime)))))/mean(Output.Shear.fullWF.TS); %Mean of the time series error of shear
    ListLidar{iNam,6} = 100*abs(Output.statistics.U.lidar.TI_mean-Output.statistics.U.fullWF.TI_mean)/Output.statistics.U.fullWF.TI_mean; % TI error
    ListLidar{iNam,7} = 100*Output.statistics.U.error.Mean_Abs_Error/Output.statistics.U.fullWF.Mean_all_full;   % mean of the error wind speed per pattern point
    ListLidar{iNam,8} = 100*abs(Output.Shear.lidar.Mean-Output.Shear.fullWF.Mean)/Output.Shear.fullWF.Mean; %  Error of the mean Shear (per simulation)
    ListLidar{iNam,9} = 100*abs(Output.REWS.lidar.mean-Output.REWS.fullWF.mean)/Output.REWS.fullWF.mean;  % Error of the mean REWS (per simulation)
    ListLidar{iNam,10} = 100*median(abs(Output.Shear.lidar.TS-Output.Shear.fullWF.TS(find(ismember(Output.TS.fullWF.time,Output.Shear.lidar.TSTime)))))/median(Output.Shear.fullWF.TS); %Mean of the time series error of shear
    ListLidar{iNam,11} = 100*mode(abs(Output.Shear.lidar.TS-Output.Shear.fullWF.TS(find(ismember(Output.TS.fullWF.time,Output.Shear.lidar.TSTime)))))/mode(Output.Shear.fullWF.TS); %Mean of the time series error of shear

%     ListLidar{iNam,8} = mean(...);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        catch
            disp(['Could not load file ' ListLidar{iNam,1} '. Continuing to the next' ])
        end
end

%% Average the seeds
ListLidarNum = cell2mat(ListLidar(:,4:end));
cnt=0;
if length(SeedVar) >1
for i = 1:length(SeedVar):size(ListLidar,1)
    cnt=cnt+1;
    LidarSeedMean{cnt,1} = ListLidar{i,1};
    LidarSeedMean{cnt,2} = ListLidar{i,2};
    LidarSeedMean{cnt,3} = ListLidar{i,3};
    ValInt = mean(ListLidarNum(i:i+length(SeedVar)-1,:));
    for ii = 1:size(ListLidarNum,2)
        LidarSeedMean{cnt,ii+3} = ValInt(ii);
    end
end
else
     LidarSeedMean = ListLidar;
end

%% Creating the contours
NumValsPlot = cell2mat(LidarSeedMean(:,2:end)); % Matrix with values, without the windfield names.

for i = 1:size(PlotNam,2)
    Var1axis = unique(NumValsPlot(:,1));
    Var2axis = unique(NumValsPlot(:,2));
    Zsurf = reshape(NumValsPlot(:,2+i),length(Var2),length(Var1));
    figure;
    contour(Var1axis,Var2axis,Zsurf,'Fill','on','Linewidth',1)
%     title(['Error REWS [%] - Vmean=' 'm/s - TI =' ' Pattern'], 'FontSize', 15)
    title(PlotNam{i},'FontSize', 15)
    xlabel(Var1Name,'FontSize', 20)
    ylabel(Var2Name,'FontSize', 20)
    grid on
    colorbar
end


