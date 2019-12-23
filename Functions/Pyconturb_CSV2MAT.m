%% Header
%
% Loading a .csv result file from PyConTurb and writting the .mat global 
% ViConDAR windfield structure file.
%
% V.Pettas/F.Costa 
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function Pyconturb_CSV2MAT(input,filename)

VariablesFile = [input.PyconturbInput_dir 'Variables_' filename '.csv'];
PyconturbFile = csvread([input.PyconturbOut_dir filename '_ConPyconturb.csv'],1); 

WF_VariablesTab = readtable(VariablesFile);
WF_Variables    = table2cell(WF_VariablesTab);
gridz = csvread([input.PyconturbInput_dir 'GridZ_' filename '.csv'],1);
gridy = csvread([input.PyconturbInput_dir 'GridY_' filename '.csv'],1);

% Parameters from PyConTurb file:
ny = length(gridy);
nz = length(gridz);
dz = abs(gridz(1,2))-abs(gridz(1,1));
dy = abs(abs(gridy(1,1))-abs(gridy(1,2)));
dt = (PyconturbFile(2,1)-PyconturbFile(1,1));
URef = WF_Variables{3};
Zh = WF_Variables{5};
gridz = (gridz)-Zh;
nt = length(PyconturbFile(:,1));
componentUinp = PyconturbFile(:,2:3:end);
componentVinp = PyconturbFile(:,3:3:end);
componentWinp = PyconturbFile(:,4:3:end);
timecomponent = PyconturbFile(:,1);

% sort matrix:
count1 = 1;
for timecount = 1:length(timecomponent)
    for col1 = 1:ny
        ComponentU{timecount}(1:ny,col1) = (componentUinp(timecount,count1:count1+(ny-1))); %#ok<*AGROW>
        ComponentV{timecount}(1:ny,col1) = (componentVinp(timecount,count1:count1+(ny-1)));
        ComponentW{timecount}(1:ny,col1) = (componentWinp(timecount,count1:count1+(ny-1)));      
        count1 = count1+ny;
    end
    ComponU(:,timecount,:) = (ComponentU{timecount})';
    ComponV(:,timecount,:) = (ComponentV{timecount})';
    ComponW(:,timecount,:) = (ComponentW{timecount})';
    count1 = 1;
end

% Grid Structure:
windfield.grid.nt = nt;
windfield.grid.ny = ny;
windfield.grid.nz = nz;
windfield.grid.dt = dt;
windfield.grid.dy = dy;
windfield.grid.dz = dz;
windfield.grid.t  = timecomponent;
windfield.grid.z  = gridz;
windfield.grid.y  = gridy;
windfield.grid.Y  = repmat(gridy,length(gridy),1);
windfield.grid.Z  = repmat(gridz',1,length(gridz));

% Windfield Structure
windfield.ny = ny;
windfield.nz = nz;
windfield.dt = dt;
windfield.u  = ComponU;
windfield.v  = ComponV;
windfield.w  = ComponW;
windfield.URef = URef;
windfield.T_offset = windfield.grid.dy*(windfield.grid.ny-1)/windfield.URef/2+windfield.grid.dt; %GridWidth/URef/2 % +dt(DS)
parsave('-v7.3',[input.PyconturbOut_dir filename '_ConPyconturb.mat'],windfield)
disp (['File ' filename '_ConPyconturb.csv succesfully converted to windfield format .mat'])
end