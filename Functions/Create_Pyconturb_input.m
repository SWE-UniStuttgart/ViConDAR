%% Header
%
% Format the output of Lidar simulator into csv readable for pyconturb
%
% V.Pettas/F.Costa 
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function Create_Pyconturb_input(Output,input,Name2Save )
% Variables:
Npoints = size(Output.Pattern.Coord,2); % number of point inth pattern
if input.nComp == 1
    Y1 = repmat(Output.Pattern.Coord(1,:),[1,3]); % temporary fix for the bug in pyconturb not accepting u only input
    Z1 = repmat(Output.Pattern.Coord(2,:)+input.Zh,[1,3]); % temporary fix for the bug in pyconturb not accepting u only input
else
    Y1 = repmat(Output.Pattern.Coord(1,:),[1,input.nComp]); % grid
    Z1 = repmat(Output.Pattern.Coord(2,:)+input.Zh,[1,input.nComp]); % grid
end
dt         = Output.Pattern.timestep_pat_vec  ; %dt of scans
TotalTime  = Output.TS.lidar.time{1, 1}  (end); % check the issue with time series (before sending the input here)
lidar_time = Output.TS.lidar.time{1, 1};

% For running PyConTurb we need the following parameters:
Variables = {TotalTime,dt,Output.statistics.U.lidar.Mean_all_LiDAR,input.nf_chunk,input.Zh,input.coh_model,input.wsp_func,input.sig_func,input.spec_func,input.seed,input.interp_data,Output.Shear.lidar.Mean,input.turb_class};
GridY     = -Output.TS.fullWF.dy*(Output.TS.fullWF.nGridY-1)/2:Output.TS.fullWF.dy:Output.TS.fullWF.dy*(Output.TS.fullWF.nGridY-1)/2;
GridZ     = (-Output.TS.fullWF.dy*(Output.TS.fullWF.nGridZ-1)/2:Output.TS.fullWF.dz:Output.TS.fullWF.dz*(Output.TS.fullWF.nGridZ-1)/2)+input.Zh;
name      = {'TotalTime','dt','u_ref','nf_chunk','z_hub','coherence_model','wsp_func','sig_func','spec_func','seed','interp_data','LiDAR_Shear','turb_class'};
Var_nam   = [name;Variables]; % creates a cell with the names and values to be written in the Variables file
warning ('off','all');
Var       = cell2table(Var_nam); 
warning ('on','all');
GridY     = array2table(GridY);
GridZ     = array2table(GridZ);

% Save Variables to .csv:
writetable(GridY,strcat(input.PyconturbInput_dir,'GridY_',Name2Save,'.csv'));
writetable(GridZ,strcat(input.PyconturbInput_dir,'GridZ_',Name2Save,'.csv'));
writetable(Var,strcat(input.PyconturbInput_dir,'Variables1_',Name2Save,'.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part removes the first line. TO DO: find a better and fster way to do it
fclose('all');
fid  = fopen( strcat(input.PyconturbInput_dir,'Variables1_',Name2Save,'.csv'),'rt');
fid2 = fopen( strcat(input.PyconturbInput_dir,'Variables_',Name2Save,'.csv'),'wt');
id = 0; %Write comments on that or remove it!
g  = fgets(fid);
while(ischar(g)) 
    id = id+1;
    if id == 1
        g = fgets(fid);
        continue
    else        
        fprintf(fid2,g);
    end
    g = fgets(fid);
end
fclose('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose number of components
if input.nComp == 1
%     compo = {'u'};
    compo = {'u','v','w'}; % temporary fix for the bug in pyconturb not accepting u only input for data_profile, dara_sig and data_spectrum options
elseif input.nComp == 2
    compo = {'u','v'};
elseif input.nComp == 3
    compo = {'u','v','w'};
end

count_ind=1;
for ind_comp=1:length(compo)
    for ind_npoints=0:(Npoints-1)
        count_ind=count_ind+1;
        columns{1,1}='index';
        columns{1,count_ind}=horzcat(compo{ind_comp},'_p',num2str(ind_npoints)); %#ok<*AGROW>
    end
end

% Creating vector with components u,v,w(k=0,1,2) required for pyconturb input con_tc
indd = 0;
if input.nComp == 1  % again  temporary fix for pyconturb bug
    for iCompVel = 0:(3-1)
        for ipoint = 1:(Npoints)
            indd = indd+1;
            k(1,indd) = iCompVel*ones(1,length(ipoint));
        end
    end
else
    for iCompVel = 0:(input.nComp-1)
        for ipoint = 1:(Npoints)
            indd = indd+1;
            k(1,indd) = iCompVel*ones(1,length(ipoint));
        end
    end
end

% Obtain matrix of time/Velocities DataFrame:
    switch input.nComp
        case 1
            for i = 1:Npoints
            VFinal_Time_U1(:,i) = Output.TS.lidar.Uval{i};
            VFinal_Time_V1(:,i) = 0.075*Output.TS.lidar.Uval{i}; % temporary fix for the bug in pyconturb not accepting u only input
            VFinal_Time_W1(:,i) = 0.05*Output.TS.lidar.Uval{i}; % temporary fix for the bug in pyconturb not accepting u only input
            end
            matrixVal = [lidar_time(1,1:end)',VFinal_Time_U1(1:end,:),VFinal_Time_V1(1:end,:),VFinal_Time_W1(1:end,:)];
        case 2
            for i = 1:Npoints
            VFinal_Time_U1(:,i) = Output.TS.lidar.Uval{i};
            VFinal_Time_V1(:,i) = Output.TS.lidar.Vval{i};
            end
            matrixVal = [lidar_time(1,1:end)',VFinal_Time_U1(1:end,:),VFinal_Time_V1(1:end,:)];
            matrixVal = matrixVal(:,any(matrixVal~=0)); % Remove columns of 0 if any
        case 3
            for i = 1:Npoints
            VFinal_Time_U1(:,i) = Output.TS.lidar.Uval{i};
            VFinal_Time_V1(:,i) = Output.TS.lidar.Vval{i};
            VFinal_Time_W1(:,i) = Output.TS.lidar.Wval{i};
            end
            matrixVal = [lidar_time(1,1:end)',VFinal_Time_U1(1:end,:),VFinal_Time_V1(1:end,:),VFinal_Time_W1(1:end,:)];
            matrixVal = matrixVal(:,any(matrixVal~=0)); % Remove columns of 0 if any
    end
MatrixVal = num2cell(matrixVal); % Create the cell required to push to csv

% Obtain matrix of spatial DataFrame:
x = (zeros(1,length(Y1))); % Pycont turb can also have shifts in X axis (parrallel to wind) (not working yet...)
for i = 1:size(columns,2)
    if i == 1
        cell1{1,i} = 'k'; % requires u=0,v=1,w=2
        cell1{2,i} = 'x';
        cell1{3,i} = 'y';
        cell1{4,i} = 'z';
    else
        jj=i-1;
        cell1{1,i} = k(jj);
        cell1{2,i} = x(jj);
        cell1{3,i} = Y1(jj);
        cell1{4,i} = Z1(jj);
    end
end

%Final cell with the  time, velocities an spatial (constraining points)
%data: Remove columns with 0
FinalCell = [columns(:,1:size(matrixVal,2));cell1(:,1:size(matrixVal,2));MatrixVal];

% Convert cell to a table and use first row as variable names
warning ('off','all');
T = cell2table(FinalCell);
warning ('on','all');
 
% Save the time constrained data frame (con_tc) to a CSV file:
writetable(T,strcat(input.PyconturbInput_dir,'con_tc1_',Name2Save,'.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part removes the first line. TO DO: find a better and faster way to do it
fclose('all');
fid  = fopen(strcat(input.PyconturbInput_dir,'con_tc1_',Name2Save,'.csv'),'rt');
fid2 = fopen(strcat(input.PyconturbInput_dir,'con_tc_',Name2Save,'.csv'),'wt');
id = 0;
a  = fgets(fid);
while(ischar(a))
    id = id+1;
    if id == 1
        a = fgets(fid);
        continue
    else       
        fprintf(fid2,a);
    end
    a = fgets(fid);
end
fclose('all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Delete useless files. This could be a problem when read/write rights are not available or when a file is open in another program
delete (strcat(input.PyconturbInput_dir,'con_tc1_',Name2Save,'.csv'))
delete (strcat(input.PyconturbInput_dir,'Variables1_',Name2Save,'.csv'))
