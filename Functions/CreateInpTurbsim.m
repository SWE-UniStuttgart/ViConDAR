%% Header
%
% Create .inp file from lidar measurements. The expressions and variables to be
% changed are hardcoipied. Then a loop goes over all the expersions and changes 
% the template file. If a new variable is added it shoyuld be added tp both expr
% and var variables.
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function CreateInpTurbsim (SaveName,TimeSeriesName,AnalysisTime,input,lidarMeas)

% List of unique expressions found in the template file 
expr={ 'UserFile','RandSeed1','NumGrid_Z','NumGrid_Y','TimeStep','AnalysisTime',...
    'UsableTime','HubHt','GridHeight','GridWidth','IECstandard','IECturbc',...
    'WindProfileType','RefHt','URef','PLExp','SCMod1','IEC_Lu',...
    'IEC_Lv','IEC_Lw'};

% List of variables to replace the expressions
var = {TimeSeriesName,num2str(input.RandSeed),num2str(lidarMeas.TS.fullWF.nGridZ),num2str(lidarMeas.TS.fullWF.nGridY),num2str(lidarMeas.Pattern.timestep_pat_vec),num2str(AnalysisTime(end)),...
    num2str(input.UsableTime),num2str(input.Zh),num2str((lidarMeas.TS.fullWF.nGridZ-1)*lidarMeas.TS.fullWF.dz),num2str((lidarMeas.TS.fullWF.nGridY-1)*lidarMeas.TS.fullWF.dy),...
    input.IECstandard,num2str(100*lidarMeas.statistics.U.lidar.TI_mean),input.WindProfileType,num2str(input.Zh),num2str(lidarMeas.statistics.U.lidar.Mean_all_LiDAR),num2str(lidarMeas.Shear.lidar.Mean),...
    input.SCMod1, input.LengthScale{1},input.LengthScale{2},input.LengthScale{3}};

% Read in the whole template and create a copy that will be changing
Template = textread('..\ConstrainedWF\Turbsim\Template_Turbsim.inp','%s','delimiter','\n'); %#ok<DTXTRD>
NewTurbSimFile = Template;

% Loop over the expressions and change one by one the values in the copy
for ii= 1:length(expr)
    expresionUsed = '"(.*?)"';
    iName  = (expr{ii});
    iVar   = var{ii};
    iLine  = find(~cellfun(@isempty,strfind(Template,(iName))));
    StrinI = Template(iLine);
    iNew   = regexprep(StrinI,expresionUsed,iVar,1);
    NewTurbSimFile{iLine} = iNew{1};
end

%Write the file
filePh = fopen(SaveName,'w');
for i = 1:size(NewTurbSimFile,1)
    fprintf(filePh,[NewTurbSimFile{i,:} '\n']);   
end
fclose(filePh);
