%% Header
%
%Dirty script to batch process full windfield statistic and create contours from
%statistics. It was done for a specific application but can be generalized
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

%%
clear all %#ok<*CLALL>
clc
addpath (genpath('..\Functions'))
% close all

%% IO definition
%Here we have to define the name of the pattern and all the expected variations

Base    = 'DTU10MW';
Shear   = 'Sh';
seed    = 'SD';
Vmean   = 'V';
TurbInt = 'TI';
PatternName = '7P_Circular';
TimestepPattern = 'Tp';
TimestepMeas   = 'Tm';

% Choose ending of file (Turbsim or Pycon)
fileEnd ='ConPyconturb_statistics';%'ConTurbSim_statistics';%

% WF statistics direcotry
WFDir = '..\ConstrainedWF\Statistics\';

%Choose which value you want to plot (index)
PlotNam = {['REWS' ' Abs Error [%] ' fileEnd] ['Shear'  ' Abs Error [%] '  fileEnd] ['Umean'  ' Abs Error [%] ' fileEnd] ['HH ' ' Abs Error [%] ' fileEnd] ['TI' ' Abs Error [%] ' fileEnd] ['U Error per slice [%] ' fileEnd] ['U Error percentage per slice [%] ' fileEnd]}; %corresponds to the outputs chosen see line 91-95

% We can combine in total 2 variations
Var1 = 5:5:25; % shear variation1
Var1Name = 'Shear';

Var2 = 5:5:30; % TI variation 2

Var3 = Var2/5;        % Tm variees only for names add the number of point in pattern here. used only when flag varPat=1
Var2Name = 'TI';
varPat =0; % flag to indicate that timestep of pattern is a variable (naming conventions)

% if you have multiple seeds you want to average:
SeedVar = 1:1;

%% Create the list of names

% if pattern timestep is a varuable usethis
if varPat== 1
    cnt=1;
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
                ListWF{cnt,1} = [WFDir curNamInt '.mat']; %#ok<*SAGROW>
                ListWF{cnt,2} = Var1(iVar1);
                ListWF{cnt,3} = Var2(iVar2);
                cnt=cnt+1;
            end
        end
    end
    
else % if you dont have timestep of pattern as a variable use this
    cnt=1;
    for iVar1 = 1:length(Var1)
        for iVar2 = 1:length(Var2)
            for iSeed = SeedVar;
                curNamInt = strrep([ Base '_' Shear num2str(Var1(iVar1),'%02.f') '_' seed num2str(SeedVar(iSeed),'%02.f') '_' Vmean '08_' TurbInt  num2str(Var2(iVar2),'%02.f') '_' PatternName '_' TimestepPattern '3d5_' TimestepMeas '00_' fileEnd],'.','d'); 
                ListWF{cnt,1} = [WFDir curNamInt '.mat'];
                ListWF{cnt,2} = Var1(iVar1);
                ListWF{cnt,3} = Var2(iVar2);
                cnt=cnt+1;
            end
        end
    end
end

%% load and save all usefull values in one cell
for iNam = 1:length(ListWF)
        try
    load(ListWF{iNam,1}); % StatisticsWF variable created // REWS Shear TI meanAbsError
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ListWF{iNam,4}  = 100*(mean(abs(StatisticsWF.Error.REWS.TS))/StatisticsWF.Original.REWS.mean);  % REWS
    ListWF{iNam,5}  = 100*(mean(abs(StatisticsWF.Error.Shear.TS))/StatisticsWF.Original.Shear.mean); % Shear
    ListWF{iNam,6}  = 100*(mean(abs(StatisticsWF.Error.Umean.TS))/StatisticsWF.Original.Umean.mean); % Umean
    ListWF{iNam,7}  = 100*(mean(abs(StatisticsWF.Error.HH.TS))/StatisticsWF.Original.HH.mean);    % HH
    ListWF{iNam,8}  = 100*abs(StatisticsWF.Error.Umean.TI)/StatisticsWF.Original.Umean.TI;       % TI
    ListWF{iNam,9}  = 100*(StatisticsWF.Error.Slice.mean)/StatisticsWF.Original.Umean.mean;       % Error per slice
    ListWF{iNam,10} = StatisticsWF.Error.Slice.mean_perc;       % Error percentage per slice

%     ListLidar{iNam,8} = mean(...);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        catch
            disp(['Could not load file ' ListWF{iNam,1} '. Continuing to the next' ])
        end
end

%% Average the seeds
ListWFNum = cell2mat(ListWF(:,4:end)); % convert to mat to apply easier numerical operations
cnt=0;
if length(SeedVar) >1
    for i = 1:length(SeedVar):length(ListWF) % loop ove seeds 
        cnt = cnt+1;      % the counter is used to produce consecutive entries to the output
        WFSeedMean{cnt,1}    = ListWF{i,1};  % save info on the case
        WFSeedMean{cnt,2}    = ListWF{i,2};
        WFSeedMean{cnt,3}    = ListWF{i,3};
        ValInt = mean(ListWFNum(i:i+length(SeedVar)-1,:)); % average the seeds
        for ii = 1:size(ListWFNum,2)
            WFSeedMean{cnt,ii+3}= ValInt(ii);  % adding to the cell the numerical valuesof each variable. Plus 3 to continue after wind field info
        end
    end
else
    WFSeedMean = ListWF;  % if no seeds just use the same table
end

%% Creating the contours
NumValsPlot=cell2mat(WFSeedMean(:,2:end)); % Matrix with permutation values and variable values, without the windfield names.

for i=1:size(PlotNam,2)
    Var1axis = unique(NumValsPlot(:,1));  % create the axes based on the unique permuatation values
    Var2axis = unique(NumValsPlot(:,2));   % create the axes based on the unique permuatation values
    % here we assume that the values are ordered according to first variable.
    % Then reshape the vector to a matrix appropriate for contouring based on the number of permutation variables
    Zsurf = reshape(NumValsPlot(:,2+i),length(Var2),length(Var1)); % here we assume that the values are ordered according to first variable. Then reshape the vector to a matrix appropriate for contouring based on the number of permutation variables
    figure;
    contour(Var1axis,Var2axis,Zsurf,'Fill','on','Linewidth',1)
%     title(['Error REWS [%] - Vmean=' 'm/s - TI =' ' Pattern'], 'FontSize', 15)
    title(PlotNam{i},'FontSize', 15,'Interpreter','None' )
    xlabel(Var1Name,'FontSize', 20)
    ylabel(Var2Name,'FontSize', 20)
    grid on
    colorbar
end


