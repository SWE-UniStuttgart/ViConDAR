%% Header
%
% Creating statistics from windfields requested. They include mean windspeed
% shear, REWS, hhTI, fullTI, hh wind speed TS. Also errors between
% the constrained and the corresponding original data are stored together. There
% is no output in the function, it just saves directly to the directory defined
% from InputParameters file.
%
% TODO: Introduce flag for calculation of frequency domain statistics to save time
%
% V.Pettas/F.Costa/M.Graefe
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function getFullWF_statistics(input,curFileInfo)

%% parameter definition
% Turbine parameters needed for REWS and shear calculations
rotor_radius = input.rotor_radius; % Radio of the Rotor [m]
Zh           = input.Zh;           % HubHeight [m]

%Processing REWS:
dist_REWS_nd = input.dist_REWS_nd; % Non dimentional span position for rotor effective wind speed calculation [define from 0 to 1 inclusive]
Wi           = input.Wi;  %Weight to be applied for rotor effective wind speed calculation
dist_REWS = rotor_radius*[dist_REWS_nd 2]; %convert ND spanwise to meters and treat the point outside with 0 weight
Wi        = [Wi 0];

%% Check if these files exist in turbsim and pycoturb folders

% get pyconturb constrained wind field names from the output folder
pyconturbWFsFull = dir([input.PyconturbOut_dir '*.mat']);
if ~isempty(pyconturbWFsFull)
    pyconturbWF_names =  extractfield(pyconturbWFsFull,'name');
    pyconWFsInd = ~cellfun(@isempty,strfind(pyconturbWF_names,curFileInfo.name)); %
    if any(pyconWFsInd)
        pyconWF = [ input.PyconturbOut_dir pyconturbWF_names{pyconWFsInd}];
        pyconWFname = pyconturbWF_names{pyconWFsInd};
    else
        pyconWF = [];
        pyconWFname = [];
    end
else
    pyconWF = [];
    pyconWFname = [];
end

% get turbsim constrained wind field names from the output folder
turbsimWFsFull = dir([input.TurbSimOut_dir '*.mat']);
if ~isempty(turbsimWFsFull)
    turbsimWF_names =  extractfield(turbsimWFsFull,'name');
    turbsimWFsInd = ~cellfun(@isempty,strfind(turbsimWF_names,curFileInfo.name)); %
    if any(turbsimWFsInd)
        turbsimWF = [input.TurbSimOut_dir turbsimWF_names{turbsimWFsInd}];
        turbsimWF_name = turbsimWF_names{turbsimWFsInd};
    else
        turbsimWF = [];
        turbsimWF_name = [];
    end
else
    turbsimWF = [];
    turbsimWF_name = [];
end

% get original wind field names from the output folder
OrWFsFull = dir([input.OriginalWF_dir '*.mat']);
if ~isempty(OrWFsFull)
    OrWF_names =  extractfield(OrWFsFull,'name');
    OrWFsInd   = ~cellfun(@isempty,strfind(OrWF_names,curFileInfo.originalWF)); %
    if any(OrWFsInd)
        ORWF = [input.OriginalWF_dir OrWF_names{OrWFsInd}];
        ORWF_name = curFileInfo.originalWF{1,1};
    else
        ORWF = [];
        ORWF_name = [];
    end
else
    ORWF = [];
    ORWF_name = [];
end

AllWf = [{ORWF} {turbsimWF} {pyconWF};{'original'} {'turbsim'} {'pyconturb'};{ORWF_name} {turbsimWF_name} {pyconWFname}];

% since we fed the function with a specific case each time, there can be maximum 3
% windfields in total here. Then we loop over the three [original
% turbsimConstrained PyconturbConstrained] wind fields and calculate statistcs
% and errors. When saving a constrained wind field the original and error
% statistics are also attached. Warning: the loop will break if there is no
% original wind field!

for iWF = 1:size(AllWf,2)
    
    %% Load windfield file
    if ~isempty(AllWf{1,iWF})
        filenameWF =AllWf{1,iWF};
        load (filenameWF); %#ok<*LOAD> % windfield variable is loaded
        
        %% Obtain and create data
        %extract components from the windfield variable
        dt       =  windfield.dt; %time step
        gridtime =  windfield.grid.nt;
        gridny   =  windfield.grid.ny;
        gridnz   =  windfield.grid.nz;
        gridz    =  -windfield.grid.z;
        gridy    =  windfield.grid.y;
        Uref     =  windfield.URef; % Mean velocity of the windfied (m/s)
        dz       =  windfield.grid.dz;
        dy       =  windfield.grid.dy;
        
        % Manipulation of data before calculations:
        compU = zeros(gridnz,gridtime,gridny);
        compV = zeros(gridnz,gridtime,gridny);
        compW = zeros(gridnz,gridtime,gridny);
        
        for i=1:gridtime
            SqueezeCompU = squeeze(windfield.u(:,i,:));
            compU(:,i,:) = flipud(SqueezeCompU');
            
            SqueezeCompV = squeeze(windfield.v(:,i,:));
            compV(:,i,:) = flipud(SqueezeCompV');
            
            SqueezeCompW = squeeze(windfield.w(:,i,:));
            compW(:,i,:) = flipud(SqueezeCompW');
        end
        
        %% Discretization
        
        fullTime       =  (dt*gridtime)-dt; %total time duration of the windfield
        fullslicesTime =  0:dt:fullTime;
        slicesDistance =  fullslicesTime*Uref; % vector with distance between slices(m)
        
        %% Calculate REWS
        
        % distances of all points in the turbsim grid from center assuming hub center =0
        for I = 1:gridny
            for II = 1:gridnz
                Distances_of_Points_In_Plane(II,I) = sqrt((gridy(I)).^2+(gridz(II)).^2); %#ok<*AGROW> %Matrix of distances from center of rotor to each point in the grid
            end
        end
        
        if input.flag_apply_weightREWS == 1
            % calculate weigths for all the grid points
            for I = 1:gridny
                for II = 1:gridnz
                    Weights_of_Points_In_Plane(II,I) = interp1(dist_REWS,Wi,Distances_of_Points_In_Plane(II,I)); %Matrix of weights for all grid points
                end
            end
            weigthTot_grid = sum(sum(Weights_of_Points_In_Plane)); %sum of weights needed for weighted average
        end
        
        %Transform Slices
        rad_values = Distances_of_Points_In_Plane<=rotor_radius;       %keep only points inside the rotor
        for ind_slicesDistance = 1:length (slicesDistance) % loop over all the slices
            % Take complete slice from turbsim:
            ExSlice_U = squeeze(compU(:,ind_slicesDistance,:)); % Values of the selected slice CompU
            ExSlice_U = ExSlice_U.*rad_values; %remove points out of rotor
            if input.flag_apply_weightREWS == 1
                ExSlice_U = ExSlice_U.*Weights_of_Points_In_Plane;         %multiply with weights
                REWS.TS(ind_slicesDistance) = sum(sum(ExSlice_U))/weigthTot_grid;
            else
                noZeroSlice = nonzeros(ExSlice_U);
                REWS.TS(ind_slicesDistance) = mean(noZeroSlice,'omitnan');
            end
        end
        REWS.mean = mean(REWS.TS,'omitnan');
        
        %% Calculate Shear power law exponent from full windfield
        %maybe add an option to calculate shear on every nth slice of the ful field to reduce time
        zero_valueZ = find(gridz==0);
        if isempty(zero_valueZ)
            [~,zero_valueZ] =min(abs(gridz));   % dirty fix for case where the grid does not include (0,0)
        end
        
        z_vec_Shear = (gridz)+Zh;       % create vector of heights
        Vhub_shear = zeros(gridtime,1); % pre assign variable
        
        % calculate power law for full turbsim data:
        for ind_sliceLaW = 1:gridtime  %for slices
            V_slice = squeeze(compU(:,ind_sliceLaW,:)); % velocity of point in slices
            Vhub_shear(ind_sliceLaW)  = mean(V_slice(zero_valueZ,:));
            v_hor   = mean (V_slice,2);     % take the average of all points in each horizontal line
            
            % find least square fit for the average vertical line
            fcn = @(alphaPL) sum((Vhub_shear(ind_sliceLaW)*(z_vec_Shear/Zh).^(alphaPL) - v_hor').^2); % least square defintion f
            [s,~,~] = fminsearch(fcn, 0.14);        % Minimise Least-Squares error
            ShearPL.TS(ind_sliceLaW) = s; %#ok<*SAGROW>
        end
        ShearPL.mean = mean(ShearPL.TS); % total shear of the wind field
        
        %% Mean wind speed
        for iTim = 1:size(compU,2)
            UmeanTs(iTim) = mean(mean(squeeze(compU(:,iTim,:))));
        end
        Umean.TS   = UmeanTs;
        Umean.mean = mean(Umean.TS);
        Umean.TI   = std(Umean.TS)/ mean(Umean.TS);
        
        %% HH wind speed TS, mean and TI
        zero_valueY = find(gridy==0); % find the 0,0 point
        zero_valueZ = find(gridz==0);
        if isempty(zero_valueZ)
            [~,zero_valueZ] = min(abs(gridz));   % dirty fix for case where the grid does not include (0,0)
        end
        if isempty(zero_valueY)
            [~,zero_valueY] = min(abs(gridy));   % dirty fix for case where the grid does not include (0,0)
        end
        
        Vhub_HH = compU(zero_valueZ,:,zero_valueY); %#ok<*FNDSB> %Velocity u at the hub height
        
        HH.TS   = Vhub_HH;
        HH.mean = mean(Vhub_HH);
        HH.TI   = std(Vhub_HH)/mean(Vhub_HH);
        
        %% Create and save .mat output
        if strcmp(AllWf{2,iWF},'original') % save only the stuff for the original if this is an orginal field
            StatisticsWF.Original.Umean = Umean;
            StatisticsWF.Original.HH    = HH;
            StatisticsWF.Original.REWS  = REWS;
            StatisticsWF.Original.Shear = ShearPL;
            
            StatisticsWF.Original.time  = fullslicesTime;
            
            StatisticsWF.Original.Parameter.rotor_radius = rotor_radius; % Radio of the Rotor [m]
            StatisticsWF.Original.Parameter.Zh = Zh;
            
            StatisticsWF.Original.grid.nGridY = gridny;
            StatisticsWF.Original.grid.nGridZ = gridnz;
            StatisticsWF.Original.grid.dy = dy;
            StatisticsWF.Original.grid.dz = dz;
            Original = StatisticsWF.Original;
            OriginalWF.U = compU;
            OriginalWF.T = fullslicesTime;
            
            %save
            if strcmp(AllWf{3,iWF}(end-3:end),'.mat')
                savename = AllWf{3,iWF}(1:end-4) ;
            else
                savename = AllWf{3,iWF};
            end
            save_data_full_path = [input.fullWF_statistics_dir savename '_statistics.mat'];
            save(save_data_full_path,'StatisticsWF')
            disp([savename ' has been processed (' datestr(datetime) '):' ])
            
            clear windfield ShearPL REWS UmeanTs Vhub_shear HH Umean StatisticsWF
        else    % save the original and the constrained data with the errors if you have a constrained windfield
            
            StatisticsWF.Constrained.Umean = Umean;
            StatisticsWF.Constrained.HH    = HH;
            StatisticsWF.Constrained.REWS  = REWS;
            StatisticsWF.Constrained.Shear = ShearPL;
            StatisticsWF.Constrained.time  = fullslicesTime;
            
            StatisticsWF.Constrained.Parameter.rotor_radius = rotor_radius;
            StatisticsWF.Constrained.Parameter.Zh = Zh;
            
            StatisticsWF.Constrained.grid.nGridY = gridny;
            StatisticsWF.Constrained.grid.nGridZ = gridnz;
            StatisticsWF.Constrained.grid.dy = dy;
            StatisticsWF.Constrained.grid.dz = dz;
            StatisticsWF.Original = Original;
            
            % Calculate errors
            Error.Umean.mean = StatisticsWF.Constrained.Umean.mean-StatisticsWF.Original.Umean.mean;
            Error.Umean.TS   = StatisticsWF.Constrained.Umean.TS -StatisticsWF.Original.Umean.TS(find(ismembertol(StatisticsWF.Original.time ,StatisticsWF.Constrained.time))); 
            Error.Umean.TI   = StatisticsWF.Constrained.Umean.TI-StatisticsWF.Original.Umean.TI;
            Error.HH.mean    = StatisticsWF.Constrained.HH.mean-StatisticsWF.Original.HH.mean;
            Error.HH.TS      = StatisticsWF.Constrained.HH.TS -StatisticsWF.Original.HH.TS(find(ismembertol(StatisticsWF.Original.time ,StatisticsWF.Constrained.time)));
            Error.HH.TI      = StatisticsWF.Constrained.HH.TI-StatisticsWF.Original.HH.TI;
            Error.Shear.mean = StatisticsWF.Constrained.Shear.mean-StatisticsWF.Original.Shear.mean;
            Error.Shear.TS   = StatisticsWF.Constrained.Shear.TS -StatisticsWF.Original.Shear.TS(find(ismembertol(StatisticsWF.Original.time ,StatisticsWF.Constrained.time)));
            Error.REWS.mean  = StatisticsWF.Constrained.REWS.mean-StatisticsWF.Original.REWS.mean;
            Error.REWS.TS    = StatisticsWF.Constrained.REWS.TS -StatisticsWF.Original.REWS.TS(find(ismembertol(StatisticsWF.Original.time ,StatisticsWF.Constrained.time)));
            % Calculate error per slice
            for iTim = 1:size(compU,2)
                [~,T_OrInd] = ismembertol(fullslicesTime(iTim),OriginalWF.T); % get the correct time matching in the original wind field
                U_Slice_ErrAll = squeeze(compU(:,iTim,:))-squeeze(OriginalWF.U(:,T_OrInd,:));
                U_Slice_ErrAllPerc = 100* U_Slice_ErrAll./squeeze(OriginalWF.U(:,T_OrInd,:));
                U_Slice_Err(iTim)  = mean(mean(abs(U_Slice_ErrAll)));
                U_Slice_ErrPerc(iTim) = mean(mean(abs(U_Slice_ErrAllPerc)));
            end
            % Frequency domain
            vWindowCon = hamming(floor(length(StatisticsWF.Constrained.time)/12)*2);
            vWindowOr  = hamming(floor(length(StatisticsWF.Original.time)/12)*2);
            dtOr       = diff(StatisticsWF.Original.time(1:2))';
            dtCon      = diff(StatisticsWF.Constrained.time(1:2))';
            
            for iY = 1:length(gridy)
                for iZ = 1:length(gridz)
                    %                     ii=ii+1;
                    %                     yCons{ii,1}=compU(iY,:,iZ);
                    %                     yOr{ii,1}=OriginalWF.U(iY,find(ismember(StatisticsWF.Original.time ,StatisticsWF.Constrained.time)),iZ);
                    [U_coh.coh(iZ,iY,:), U_coh.f(iZ,iY,:)]        = mscohere(OriginalWF.U(iZ,find(ismembertol(StatisticsWF.Original.time ,StatisticsWF.Constrained.time)),iY),compU(iZ,:,iY),round(size(compU,2)/20),[],[],1/(StatisticsWF.Constrained.time(2)-StatisticsWF.Constrained.time(1)));
                    [Fdom.Or.S(iZ,iY,:),Fdom.Or.f(iZ,iY,:)]       = pwelch(detrend(OriginalWF.U(iZ,:,iY),'constant'),vWindowOr,[],[],1/dtOr,'onesided');
                    [Fdom.Const.S(iZ,iY,:),Fdom.Const.f(iZ,iY,:)] = pwelch(detrend(compU(iZ,:,iY),'constant'),vWindowCon,[],[],1/dtCon,'onesided');
                end
            end
            %             warning('off','all')
            %             [Sxy,fYiyin] = CalcEstimateCPSD(yCons,yOr,'fs',0.5,'nFFT',202);
            %             [SCon,~] = CalcEstimateCPSD(yCons,yCons,'fs',0.5,'nFFT',202);
            %             [SOr,~] =CalcEstimateCPSD(yOr,yOr,'fs',0.5,'nFFT',202);
            %             warning('on','all')
            %             U_coh.YYmean_cxy = abs(Sxy).^2./abs(SCon)./abs(SOr);
            %             U_coh.YYFmean =fYiyin;
            %             [r,lags] =xcorr(OriginalWF.U(iY,:,iZ),compU(iY,:,iZ),20);
            %             [rOr,lagsOr] =xcorr(OriginalWF.U(iY,find(ismember(StatisticsWF.Original.time ,StatisticsWF.Constrained.time)),iZ),20);
            %             [rCon,lagsCon] =xcorr(compU(iY,:,iZ),20);
            U_coh.CohMean = squeeze(mean(mean(U_coh.coh,1)) );
            U_coh.fMean   = squeeze(U_coh.f(1,1,:));
            Error.Slice.TS   = U_Slice_Err;
            Error.Slice.mean = mean(U_Slice_Err);
            Error.Slice.TS_perc   = U_Slice_ErrPerc;
            Error.Slice.mean_perc = mean(U_Slice_ErrPerc);
            StatisticsWF.Error = Error;
            StatisticsWF.U_coh = U_coh;
            StatisticsWF.Fdom  = Fdom;
            %save
            if strcmp(AllWf{3,iWF}(end-3:end),'.mat')
                savename = AllWf{3,iWF}(1:end-4) ;
            else
                savename = AllWf{3,iWF};
            end
            save_data_full_path = [input.fullWF_statistics_dir savename '_statistics.mat'];
            save(save_data_full_path,'StatisticsWF')
            disp([savename ' has been processed (' datestr(datetime) '):' ])
        end
    else
        disp(['No file with name ' curFileInfo.name ' found in ' AllWf{2,iWF} ' folder'])
        
    end
    clear windfield ShearPL REWS UmeanTs Vhub_shear HH Umean StatisticsWF
end
