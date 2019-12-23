%% Header
%
% Dirty script to batch process full windfield statistics and compare time series
% betweeen constrained and original windfields. It can handel one original wind
% field with may variations ofthe constraining (eg Pyconturb,turbsim, patterns, Tp and Tm)
% It was done for a specific application but can be generalized
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

clearvars
clc
close all
addpath (genpath('..\Functions'))

% direc=fullfile('X:','ViConDAR_Test_Hor_PL','ConstrainedWF','Statistics\');
direc = '..\ConstrainedWF\Statistics\';
filesAll = dir(fullfile(direc, '*.mat'));
files    = extractfield(filesAll,'name')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define names as in inputs!
input.nameBase = 'DTU10MW';

% Don't use more than one original wind field here. One value per variable only!
input.freeInp = {'Sh' ,[5];
    'SD' ,[1];
    'V'  ,[8];
    'TI' ,[5];
    }; %#ok<*NBRAK>

% The lidar parameters can vary freely!!!
% parameters of the simulator cannot be changed
input.fixedInp={
    'Pat' ;   %input.PatternNames
    %                 'Ns';     %input.noise_U
    'Tp';     %timeStep_pattern
    'Tm'   ;  %input.timeStep_Measurements
    %                 'Pos'   ; %input.Pos_LiDAR
    %                 'Fd'   ;  %input.ref_plane_dist
    %                 'DAv'  ;  %input.distance_av_space
    %                 'SlAv'  ; %points_av_slice
    };
input.PatternNames = {'7P_Circular' };  % names of the patterns. Important: number of names should equal number of Y,Z coordinates
input.timestep_pat_vec      = {[2 ]  }; %Time step of the total pattern. Sampling rate of total pattern should be that npoins*timestep_meas<=timestep_pat(s). Add one value for each pattern
input.timeStep_Measurements = {[0] }; %Time step between each single measured point. Add one value for each pattern [s]
input.ref_plane_dist = [250];   % Reference Plane for LOS (distance[m])
input.Pos_LiDAR      = [0,0]; % LiDAR position offsetfrom hub center(meters)==> [Y,Z] WE NEED TO FIX THE APPLICATION OF OFFSET IN LOS ONLY!!!! It cannot be used to loop over it. It has to be fixed for now
input.distance_av_space = [0];    % [m] values to use for imitating range gate averaging in the calcualtion of wind speeds from pulses meters ahead and afer the point
input.points_av_slice   = [0];     % how many point/slices you want to take in the averaging of distance_av_slice  Totalpoints = distance_av_slice/points_av_slice+1 IT HAS TO BE AN EXACT DIVISION FOR NOW!!!!
input.noise_U  = [20];
input.AllFixed = {'Pat';'Ns';'Tp';'Tm';'Pos';'Fd' ;'DAv';'SlAv'}; % This should not be changed. Only when a new lidar feature is added. This is a list of all the parameters requires to run lidar simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

perm_cell = getNamesFromInputs(input);
NamesInit = perm_cell.OutNames;

Errcell = {'Error REWS' 'Error Umean' 'Error Shear'};
VarCell = {'Umean' 'Shear' 'HH'};
% We can combine in total 2 variations
Var1 = 15; % shear variation1
Var1Name = 'Shear';

Var2 = 5; % TI variation 2
Var3 = Var2/5;        % Tm variees only for names add the number of point in pattern here
Var2Name = 'TI';
varPat = 0; % flag to indicate that timestep of pattern is a variable (naming conventions)

% if you have multiple seeds you want to average:
SeedVar = 1;
%%
cunt=0;
for nWF = 1:size(NamesInit,1) % for each WF it arranges in different .structures for each time step
    
    curNam   = NamesInit{nWF,1};
    filesInd = ~cellfun(@isempty,strfind(files,curNam)); %
    FileNam  = files(filesInd);
    for iNam = 1:size(FileNam)
        cunt = cunt+1;
        load([direc FileNam{iNam}]); % Load each loop one statistic
        WFdata{cunt,1} = FileNam{iNam};
        WFdata{cunt,2} = perm_cell.values(nWF);
        WFdata{cunt,3} = []; % Umean constrained
        WFdata{cunt,4} = StatisticsWF.Constrained.Umean.TS     ; % Umean constrained
        WFdata{cunt,5} = StatisticsWF.Original.Umean.TS     ; % Umean original
        WFdata{cunt,6} = StatisticsWF.Constrained.Shear.TS    ; %Shear constrained
        WFdata{cunt,7} = StatisticsWF.Original.Shear.TS    ;  %Shear original
        WFdata{cunt,8} = StatisticsWF.Constrained.HH.TS; % %HH constrained
        WFdata{cunt,9} = StatisticsWF.Original.HH.TS; % %HH original
        WFdata{cunt,10} = abs(StatisticsWF.Error.REWS.TS); % Error  REWS
        WFdata{cunt,11} = abs(StatisticsWF.Error.Umean.TS); % %Error Umean StatisticsWF.Error.Shear.TS
        WFdata{cunt,12} = abs(StatisticsWF.Error.Shear.TS); % %Error UShear
        WFdata{cunt,13} = StatisticsWF.Error.Slice.TS; % %Error UShear  
        WFdata{cunt,14} = abs(StatisticsWF.Error.Shear.TS); % %Error UShear  
        WFdata{cunt,15} = StatisticsWF.Error.Slice.TS_perc; % %Error Slice       
        WFdata{cunt,16} = StatisticsWF.Original.time; % %Time original
        WFdata{cunt,17} = StatisticsWF.Error.Slice.TS_perc; % %Time constrained
    end
end

%% PLOTTING

for i = 10:size(WFdata,2)-2 % loop over variables
    legCell = {};
    %plot Umean error
    figure
    for iCase2 = 1:size(WFdata,1) %loop over different cases for the group of cases
        curNam = WFdata{iCase2};
        plot(WFdata{iCase2,17},WFdata{iCase2,i},'LineWidth',1.5)
        legCell{1,end+1} = curNam; %#ok<*SAGROW>
        grid on
        hold on
        
    end
    xlabel('Time [s]')
    ylabel(Errcell{i-9})
    set(gca,'FontSize',14)
    legend(legCell,'Interpreter','None','FontSize',10)
    clear curNam
    hold off
end

cunt3 = 0;
for i = 4:2:9 % loop over variables
    cunt3 = cunt3+1;
    legCell2={};
    %plot Umean error
    figure
    for iCase2 = 1:size(WFdata,1) %loop over different cases for the group of cases        
        curNam = WFdata{iCase2};      
        plot(WFdata{iCase2,17},WFdata{iCase2,i},'LineWidth',1.5)
        legCell2{1,end+1} = curNam;
        grid on
        hold on       
    end
    plot(WFdata{iCase2,16},WFdata{iCase2,i+1},'LineWidth',1.5)
    legCell2{1,end+1} = 'Original';
    xlabel('Time [s]')
    ylabel(VarCell{cunt3})
    set(gca,'FontSize',14)
    legend(legCell2,'Interpreter','None','FontSize',10)
    clear curNam
    hold off
end
