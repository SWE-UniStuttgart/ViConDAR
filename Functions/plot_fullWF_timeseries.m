%% Header
%
% Plotting the time series from full wind fields at grid points requested by the
% user
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function plot_fullWF_timeseries(input,perm_cell)

pyconturbWFsFull = dir([input.PyconturbOut_dir '*.mat']);
if ~isempty(pyconturbWFsFull)
    pyconturbWF_names =  extractfield(pyconturbWFsFull,'name');
end

turbsimWFsFull = dir([input.TurbSimOut_dir '*.mat']);
if ~isempty(turbsimWFsFull)
    turbsimWF_names =  extractfield(turbsimWFsFull,'name');
end

for iPoint = 1:size(input.points_plot_WF_timeseries,2) % loop over requested points in space
    Ypoint = input.points_plot_WF_timeseries(1,iPoint); %Y point in meters
    Zpoint = input.points_plot_WF_timeseries(2,iPoint); %Z  point in meters    
    
    %Go through orginal wind field folders to match the names
    for iNamOr = 1: length(perm_cell.namesOWF) % loop over orifinal windfields
        figure
        legCell = {};
        curNamOr = perm_cell.namesOWF{iNamOr,1};
        load([input.OriginalWF_dir curNamOr]); % variable windfield loaded
        % find the closesrt point in the grid as requeted by the user
        YpointGrid = find(abs(Ypoint-windfield.grid.y) == min(abs(Ypoint-windfield.grid.y)),1); %Y point as index
        ZpointGrid = find(abs(Zpoint-windfield.grid.z) == min(abs(Zpoint-windfield.grid.z)),1); %z point as index
        plot(windfield.grid.t,windfield.u(YpointGrid,:,ZpointGrid),'LineWidth',1.5)
        legCell{1,end+1} = curNamOr; %#ok<*AGROW>
        hold on
        
        % loop over constrained wind fields with pyconturb that come from the same
        % original windfield
        if input.plot_WF_timeseries_Pyconturb == 1 ;
            pyconWFsInd = ~cellfun(@isempty,strfind(pyconturbWF_names,curNamOr)); %
            pyconWFs = pyconturbWF_names(pyconWFsInd);
            for ipyconWF = 1:size(pyconWFs,2)
                curNamPycon = pyconWFs{ipyconWF};
                load([input.PyconturbOut_dir curNamPycon]); % variable windfield loaded
                % find the closesrt point in the grid as requeted by the user
                YpointGrid = find(abs(Ypoint-windfield.grid.y)==min(abs(Ypoint-windfield.grid.y)),1); %Y point as index
                ZpointGrid = find(abs(Zpoint-windfield.grid.z)==min(abs(Zpoint-windfield.grid.z)),1); %z point as index
                plot(windfield.grid.t,windfield.u(YpointGrid,:,ZpointGrid),'LineWidth',1.5)
                legCell{1,end+1}= curNamPycon(1:end-4);
                hold on
            end
        end
        
        % loop over constrained wind fields with turbsim that come from the same
        % original windfield
        if input.plot_WF_timeseries_Turbsim == 1 ;
            turbsimWFsInd = ~cellfun(@isempty,strfind(turbsimWF_names,curNamOr)); %
            turbsimWFs = turbsimWF_names(turbsimWFsInd);
            for iturbsimWF = 1:size(turbsimWFs,2)
                curNamTurbsim = turbsimWFs{iturbsimWF};
                load([input.TurbSimOut_dir curNamTurbsim]); % variable windfield loaded
                % find the closesrt point in the grid as requeted by the user
                YpointGrid = find(abs(Ypoint-windfield.grid.y) == min(abs(Ypoint-windfield.grid.y)),1); %Y point as index
                ZpointGrid = find(abs(Zpoint-windfield.grid.z) == min(abs(Zpoint-windfield.grid.z)),1); %z point as index
                plot (windfield.grid.t,windfield.u(YpointGrid,:,ZpointGrid),'LineWidth',1.5)
                legCell{1,end+1} = curNamTurbsim(1:end-4);
                hold on
            end
        end    
        grid on
        title (['Timeseries at Y:' num2str(Ypoint) 'm Z:' num2str(Zpoint) 'm'])
        xlabel ('Time [s]')
        ylabel ('U vel [m/s]')
        set(gca,'FontSize',14)
        legend (legCell,'Interpreter','None','FontSize',10)
        clear legCell
    end
end
