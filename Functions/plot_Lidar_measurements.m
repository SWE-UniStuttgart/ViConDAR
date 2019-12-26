%% Header
%
% Plotting the lidar outputs. Comparison of original time series and lidar
% measurements and convergence metrics.
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function plot_Lidar_measurements(input,perm_cell)

% get names of files with lidar outputs that have the same name base, free inputs
% and pattern to plot them together
PatInd = find(strcmp((perm_cell.variables{1, 1}),'Pat')==1);
FreeInpInd = size(input.freeInp,1);
MatchInd = [1:FreeInpInd PatInd];

%find all combination of same inputs
for iRow1 = 1:size(perm_cell.values,1)
    for iRow2 = 1:size(perm_cell.values,1)
        Comb(iRow1,iRow2) = isequal(cell2mat(perm_cell.values {iRow1,1}(1,MatchInd)),cell2mat(perm_cell.values{iRow2,1}(1,MatchInd)));
    end
end
TotalComb = unique(Comb,'rows'); % get only the unique combinations

% Cases that should be plot together:
for iComb = 1:size(TotalComb,1)
    iTotalComb  = (TotalComb(iComb,:));
    Case{iComb} = find(iTotalComb==1); %#ok<*AGROW>
end

for iCase = 1:size(Case,2) %loop over the groups of cases
    curCase = Case{iCase};
    % plot all lidar points
    for iPoint = 1:length(input.PatternY{perm_cell.values{curCase(1),1}{1,PatInd},1}) % loop over pattern points
        figure
        for iCase1 = 1:length(curCase) %loop over different cases for the same point in the pattern
            curNam = perm_cell.OutNames{curCase(iCase1),1};
            try
                load ([input.LidarOutput_dir curNam '.mat'])
            catch
                error(['Could not find lidar ouptut file ' input.LidarOutput_dir curNam '.mat'])
            end
            if iCase1 == 1
                %  plot(Output.TS.fullWF.time,Output.TS.fullWF.Uval{1, ii})
                plot(Output.TS.fullWF.time,Output.TS.fullWF.Uval{1, iPoint}, Output.TS.lidar.time{1, iPoint},Output.TS.lidar.Uval{1, iPoint},'LineWidth',1.5)
                legCell = {'Original WF', curNam};
            else
                plot(Output.TS.lidar.time{1, iPoint},Output.TS.lidar.Uval{1, iPoint},'LineWidth',1.5)
                legCell{1,end+1} = curNam;
            end
            grid on
            hold on
        end
        title (['Pattern' input.PatternNames{perm_cell.values{curCase(1),1}{1,PatInd}} ' Y:' num2str(input.PatternY{perm_cell.values{curCase(1),1}{1,PatInd},1}(iPoint)) 'm Z:' num2str(input.PatternZ{perm_cell.values{curCase(1),1}{1,PatInd},1}(iPoint)) 'm'], 'Interpreter', 'none')
        xlabel ('Time [s]')
        ylabel ('U vel [m/s]')
        set(gca,'FontSize',14)
        legend (legCell,'Interpreter','None','FontSize',10)
        hold off
        clear legCell
    end
    clear curNam
    
    % plot Shear
    figure
    for iCase2 = 1:length(curCase) %loop over different cases for the group of cases
        curNam = perm_cell.OutNames{curCase(iCase2),1};
        load ([input.LidarOutput_dir curNam '.mat'])
        if iCase2 == 1
            plot(Output.TS.fullWF.time,Output.Shear.fullWF.TS,Output.Shear.lidar.TSTime,Output.Shear.lidar.TS,'LineWidth',1.5)
            legCell = {'Original WF', curNam};
        else
            plot(Output.Shear.lidar.TSTime,Output.Shear.lidar.TS,'LineWidth',1.5)
            legCell{1,end+1} = curNam;
        end      
        grid on
        hold on
    end
    title (['Pattern' input.PatternNames{perm_cell.values{curCase(1),1}{1,PatInd}}],'Interpreter','None')
    xlabel ('Time [s]')
    ylabel ('Shear exponent [-]')
    set(gca,'FontSize',14)
    legend (legCell,'Interpreter','None','FontSize',10)
    clear curNam
    hold off
    
    % plot REWS
    figure
    for iCase2 = 1:length(curCase) %loop over different cases for the group of cases
        curNam = perm_cell.OutNames{curCase(iCase2),1};
        load ([input.LidarOutput_dir curNam '.mat'])
        if iCase2 == 1
            plot(Output.TS.fullWF.time,Output.REWS.fullWF.TS,Output.REWS.lidar.TSTime,Output.REWS.lidar.TS,'LineWidth',1.5)
            legCell = {'Original WF', curNam};
        else
            plot(Output.REWS.lidar.TSTime,Output.REWS.lidar.TS,'LineWidth',1.5)
            legCell{1,end+1} = curNam;
        end       
        grid on
        hold on
    end
    title (['Pattern' input.PatternNames{perm_cell.values{curCase(1),1}{1,PatInd}}],'Interpreter','None')
    xlabel ('Time [s]')
    ylabel ('REWS [m/s]')
    set(gca,'FontSize',14)
    legend (legCell,'Interpreter','None','FontSize',10)
    clear curNam curCase
    hold off
end
