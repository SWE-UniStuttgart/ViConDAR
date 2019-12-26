%% Header
%
% WRAPPER script for ViConDAR.
% Processing and Postprocessing data from Lidar Simulator, comparison and
% constraining Windfields. All the paramaters and flags are assigned in the
% InputParameters function.
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

clc;
close all;
clearvars; 
tmp = matlab.desktop.editor.getActive;
cd(fileparts(tmp.Filename));
% Add functions, Main and HelpfulStandAlone folders to the path
addpath(genpath('../Functions'))
addpath(genpath('../Main'))
addpath(genpath('../HelpfulStandAlone'))

%% Input variables and directory creation
input = InputParameters(); % Get all the inputs
perm_cell = getNamesFromInputs(input); % Reading the InputParameters and assign parameters for each name permutation

% Creating directories if they don't exist (here because input.<path> should exist):
for indDirec = 1:length(input.AddUserDir)
    FolderName = input.AddUserDir{indDirec};
    if ~exist(FolderName,'dir')
        mkdir(FolderName);
    end
end

%% Obtain the lidar output form virtual lidar
if input.flag_getLidarOutput == 1
    % Loop over requested wind fields
    for iProcess = 1:size(perm_cell.OutNames,1)
        % pass information to getLidarOutput
        curFileInfo.name   = perm_cell.OutNames{iProcess};
        curFileInfo.values = perm_cell.values{iProcess};
        curFileInfo.variables = perm_cell.variables;
        
        % find the corresponding original WF based on names
        for iName = 1:size(perm_cell.namesOWF,1)
            Indexi = strfind (perm_cell.OutNames{iProcess},perm_cell.namesOWF{iName});
            if Indexi == 1
                indexvec(iName) = 1;
            else
                indexvec(iName) = 0;
            end
        end
        curFileInfo.originalWF = perm_cell.namesOWF (find(indexvec==1))  ; %#ok<*FNDSB>
        Output = getLidarOutput(input,curFileInfo); % Obtain lidar measurements
        disp([curFileInfo.name ' has been processed (' datestr(datetime) '):' ])
    end
    disp('Creating virtual lidar output finished successfully')
end

%% Obtain the inputs for constraining WF with turbsim

if input.flag_getTurbsimInput == 1
    for iProcess = 1:size(perm_cell.OutNames,1)
        % pass information to turbsim input
        curFileInfo.name   = perm_cell.OutNames{iProcess};
        curFileInfo.values = perm_cell.values{iProcess};
        curFileInfo.variables = perm_cell.variables  ;
        try
            load([input.LidarOutput_dir curFileInfo.name]); % Output file is loaded here!
        catch e
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'There was an error! The message was:\n%s',e.message);
            error(['Error loading lidar measurements for creating Turbsim inputfile ' curFileInfo.name '. Check if the lidar measurement exists and is in the correct folder.'])
        end
        [VFinal,t_lidar,Y1,Z1] = Create_Turbsim_Vinput(Output,input); % Creating time series and information needed for turbsim timeSR file 
        [~,RefNode] = min(abs(Z1-input.Zh)+abs(Y1));
        writeTurbSimTimeSeriesInput([input.TurbSimInput_dir curFileInfo.name '.TimeSer'], curFileInfo.name, VFinal, t_lidar', Y1,Z1,input.nComp,RefNode); % Write time series file for TurbSim
        CreateInpTurbsim([input.TurbSimInput_dir curFileInfo.name '_ConTurbSim.inp'],[curFileInfo.name '.TimeSer'],t_lidar,input,Output); %create .inp file for turbsim based on the template
        disp([curFileInfo.name ' input for turbsim done.' ])
    end
end

%% Obtain the inputs for constraining WF with pyconturb

if input.flag_getPyconturbInput==1
    for iProcess = 1:size(perm_cell.OutNames,1) % Get data from all requested wind fields. Based on names it takes variables and values.
        curFileInfo.name   = perm_cell.OutNames{iProcess}  ;
        curFileInfo.values = perm_cell.values{iProcess}  ;
        curFileInfo.variables = perm_cell.variables  ;
        disp('Starting creating input files for PyConTurb')
        try
            load([input.LidarOutput_dir curFileInfo.name]); % Output file is loaded here!
        catch e
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'There was an error! The message was:\n%s',e.message);
            error (['Error loading lidar measurements for creating Pyconturb input file ' curFileInfo.name '. Check if the lidar measurement exists and is in the correct folder.'])
        end
        Create_Pyconturb_input(Output,input,curFileInfo.name)% writing the csv files con_tc,variables,gridY and gridZ we need to run pyconturb
        disp([curFileInfo.name ' input for PyConTurb done'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate automatically resample factor:Wr
%         format long
%         Vector_resampling_factor=2:1:8;% Vector of possible resampling factor values
%         original_TimeStep=.125;% time step of the original windfield
%         for prueb_fact_ind=1:length(Vector_resampling_factor);
%             number_fact(prueb_fact_ind)=timestep_pat/Vector_resampling_factor(prueb_fact_ind); %#ok<*SAGROW>
%             number_fact(prueb_fact_ind)=round(number_fact(prueb_fact_ind),3);
%             if number_fact(prueb_fact_ind)<original_TimeStep
%                 break  %%% Esta parte hay que mejorarla.  No usar break!1
%             else
%                 busquIndex=TOTAL_TIME_SERIE_DURATION/ number_fact(prueb_fact_ind);
%                 mod(busquIndex,1);
%                 if mod(busquIndex,1)==0
%                     RF(prueb_fact_ind)=Vector_resampling_factor(prueb_fact_ind);
%                 else
%                     RF(prueb_fact_ind)=0;
%                 end
%                 resampling_factor=max(RF);
%                 index_resamplingFactor=find(Vector_resampling_factor==resampling_factor);
%                 value_newSteptime=number_fact(index_resamplingFactor);
%             end
%         end
%         format short
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Obtain constrained windfield with Turbsim

if input.flag_obtain_Con_Turbsim == 1
    for iProcess = 1:size(perm_cell.OutNames,1)  % try parfor here 
        curFileInfo.name   = perm_cell.OutNames{iProcess}; 
        curFileInfo.values = perm_cell.values{iProcess};
        curFileInfo.variables = perm_cell.variables;
        filenamTurbSimConInp  =  [curFileInfo.name]; 
        RunTurbSimCon_WF(input,filenamTurbSimConInp) % runs turbsim with system command
    end
end

%% Obtain constrained windfield with Pyconturb running from matlab

if input.flag_obtain_Con_Pyconturb_Matlab ==1
    for iProcess = 1:size(perm_cell.OutNames,1)
        curFileInfo.name = perm_cell.OutNames{iProcess};
        RunPyconturb_Con_WF(input,curFileInfo.name) % requires python.exe path with PyConturb compatible environment
    end
end

%% Obtain constrained windfield with Pyconturb running externally with python wrapper

if input.flag_obtain_Con_Pyconturb_python == 1
    for iProcess = 1:size(perm_cell.OutNames,1)
        curFileInfo.name = perm_cell.OutNames{iProcess};
        %create names for the python wrapper
        Pythonwrapperin{iProcess,1} = [input.PyconturbInput_dir 'con_tc_' curFileInfo.name '.csv,' input.PyconturbInput_dir 'GridY_' curFileInfo.name '.csv,' input.PyconturbInput_dir 'GridZ_' curFileInfo.name '.csv,' input.PyconturbInput_dir 'Variables_' curFileInfo.name '.csv,' input.PyconturbOut_dir curFileInfo.name '_ConPyconturb.csv']; %#ok<*SAGROW>
    end
    % write a .csv file with all the inputs for the python wrapper. It will
    % be read from PyConTurb. Wrapper.csv is saved in the same folder as the
    % rest of .csv files
    disp('Starting creating input for python wrapper for PyConTurb')
    filePh = fopen([input.PyconturbInput_dir 'PythonWrapperInput.csv'],'w');
    for iNamCsv = 1:size(Pythonwrapperin,1)+1
        if   iNamCsv == 1
            fprintf(filePh,'%s\n','contc,gridy,gridz,variables,savename');
        else
            fprintf(filePh,'%s\n',Pythonwrapperin{iNamCsv-1,1});
        end
    end
    fclose(filePh);
    disp('Input for python wrapper for PyConTurb created succesfully')
    % After having finished the simulation in python this part translates CSV
    % to global windfield .mat format
    if input.flag_obtain_Con_Pyconturb_ConverToMat == 1
        for iProcess = 1:size(perm_cell.OutNames,1)
            curFileInfo.name = perm_cell.OutNames{iProcess};
            Pyconturb_CSV2MAT(input,curFileInfo.name)
        end
    end
end

%% Calculate statistics from full wind fields
if input.flag_calculate_fullWF_statistics == 1
    disp('Starting calculations of statistics for full wind fields ')
    for iProcess = 1:size(perm_cell.OutNames,1)
        % pass required information about windfields
        curFileInfo.name   = perm_cell.OutNames{iProcess}  ;
        curFileInfo.values = perm_cell.values{iProcess}  ;
        curFileInfo.variables = perm_cell.variables  ;
        
        % find the corresponding original WF 
        for iName = 1:size(perm_cell.namesOWF,1)
            Indexi = strfind (perm_cell.OutNames{iProcess},perm_cell.namesOWF{iName});
            if Indexi == 1
                indexvec(iName) = 1;
            else
                indexvec(iName) = 0;
            end
        end
        curFileInfo.originalWF = perm_cell.namesOWF (find(indexvec==1))  ;
        getFullWF_statistics(input,curFileInfo); % Reading the file names and extracting statistics from the result files(global .mat windfields)
        clear indexvec
    end
    disp('Calculations of statistics for full wind fields done')
end
    
%% Plotting Lidar measurements vs original windfield
% the requested cases will be plotted
if input.flag_plot_lidar == 1
    disp('Starting plotting lidar measurements and statistics')
    plot_Lidar_measurements(input,perm_cell)
    disp('Plotting lidar measurements and statistics done')
end

%% Plotting time series from windfield (e.g. constrained vs original)
if input.flag_plot_WF_timeseries  == 1
    disp('Starting plotting requested full wind field timeseries')
    plot_fullWF_timeseries(input,perm_cell)
    disp('Plotting requested full wind field timeseries done')
end

%% Plotting contours from windfield slices  (e.g. constrained vs original)
if input.plot_fullWF_Slices == 1
    disp('Starting plotting requested full wind field slices')
    plot_fullWF_Slices(input,perm_cell)
    disp('Plotting requested full wind field slices done')
end
