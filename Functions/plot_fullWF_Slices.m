%% Header
%
% Plotting requested slices from full wind fields at timesteps requested by the
% user as defined in input.time_fullWF_Slices.
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function plot_fullWF_Slices(input,perm_cell)

disp('Starting plotting full wind field slices requested')
pyconturbWFsFull = dir([input.PyconturbOut_dir '*.mat']);
if ~isempty(pyconturbWFsFull)
    pyconturbWF_names =  extractfield(pyconturbWFsFull,'name');
end

turbsimWFsFull = dir([input.TurbSimOut_dir '*.mat']);
if ~isempty(turbsimWFsFull)
    turbsimWF_names = extractfield(turbsimWFsFull,'name');
end

for iTimeSlice = 1:length(input.time_fullWF_Slices) % loop over requested slices in time
    curTime = input.time_fullWF_Slices(iTimeSlice); % current requested time
    
    %Go through original windfield folder to match the names 
    for iNamOr = 1: length(perm_cell.namesOWF) % loop over orifinal windfields
        figure
        curNamOr = perm_cell.namesOWF{iNamOr,1};
        load([input.OriginalWF_dir curNamOr]); % variable windfield loaded
        % find the closest timestep in the grid as requeted by the user
        TSlice = find(abs(curTime-windfield.grid.t)==min(abs(curTime-windfield.grid.t)),1); %Y point as index
        contour (windfield.grid.y,windfield.grid.z,(squeeze(windfield.u(:,1,:)))','Fill','on')
        colorbar('peer',gca);
        CLimits = [floor(min(min(squeeze(windfield.u(:,TSlice,:))))) ceil(max(max(squeeze(windfield.u(:,TSlice,:)))))];
        set(gca,'BoxStyle','full','CLim',CLimits,'Layer','top');
        xlabel ('Y axis [m]')
        ylabel ('Z axis [m]')
        set(gca,'FontSize',14)
        title (['U comp at t=:' num2str(windfield.grid.t(TSlice)) ' s WF:' curNamOr],'Interpreter','None','Interpreter','None','FontSize',10)
        
        % loop over constrained wind fields with pyconturb that come from the same
        % original windfield
        if input.plot_WF_slices_Pyconturb == 1 ;
            pyconWFsInd = ~cellfun(@isempty,strfind(pyconturbWF_names,curNamOr)); %
            pyconWFs    = pyconturbWF_names(pyconWFsInd);
            for ipyconWF = 1:size(pyconWFs,2)
                curNamPycon = pyconWFs{ipyconWF};
                load([input.PyconturbOut_dir curNamPycon]); % variable windfield loaded
                % find the closesrt timestep in the grid as requeted by the user
                TSlice = find(abs(curTime-windfield.grid.t) == min(abs(curTime-windfield.grid.t)),1); %Y point as index
                figure
                contour(windfield.grid.y,windfield.grid.z,(squeeze(windfield.u(:,TSlice,:)))','Fill','on')
                colorbar('peer',gca);
                xlabel ('Y axis [m]')
                ylabel ('Z axis [m]')
                set(gca,'FontSize',14)
                title (['U comp at t=' num2str(windfield.grid.t(TSlice)) ' s WF:' curNamPycon],'Interpreter','None','FontSize',10)
            end
        end
        
        % loop over constrained wind fields with turbsim that come from the same
        % original windfield
        if input.plot_WF_slices_Turbsim == 1 ;
            turbsimWFsInd = ~cellfun(@isempty,strfind(turbsimWF_names,curNamOr)); %
            turbsimWFs = turbsimWF_names(turbsimWFsInd);
            for iturbsimWF = 1:size(turbsimWFs,2)
                curNamTurbsim = turbsimWFs{iturbsimWF};
                load([input.TurbSimOut_dir curNamTurbsim]); % variable windfield loaded
                % find the closesrt timestep in the grid as requeted by the user
                TSlice = find(abs(curTime-windfield.grid.t) == min(abs(curTime-windfield.grid.t)),1); %Y point as index
                figure
                contour (windfield.grid.y,windfield.grid.z,(squeeze(windfield.u(:,TSlice,:)))','Fill','on')
                colorbar('peer',gca);
                set(gca,'BoxStyle','full','CLim',CLimits,'Layer','top');
                xlabel ('Y axis [m]')
                ylabel ('Z axis [m]')
                set(gca,'FontSize',14)
                title (['U comp at t=' num2str(windfield.grid.t(TSlice)) ' s WF:' curNamTurbsim],'Interpreter','None','FontSize',10)
                
            end
        end
    end
    clear CLimits
end