%% Header
%
% Once we have obtained the inputt for TurbSim, this function
% runs turbsim to obtain the constrained windfield in .wnd format. After the 
% simulation is done it converts the .wnd to global .mat format.
%
% V.Pettas/F.Costa 
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

function RunTurbSimCon_WF(input,Filename)

inputfile = [Filename '_ConTurbSim' '.inp'];

% Copy required files where Turbsim since it cannot work with relative paths
destination = input.pathToTurbSim  ;
origen1 = [input.TurbSimInput_dir inputfile]; % copying input file to turbsim folder
origen2 = [input.TurbSimInput_dir Filename '.TimeSer']; % copying timeseries file to turbsim folder
try
    copyfile(origen1,destination);
    copyfile(origen2,destination);
catch
    error(['Error loading Turbsim inputs ' Filename ' for contrained turbulence. Check whether the inputs exist and are in the correct folder.'])
end

[status,cmdout] = system([input.TurbSimExe_path ' ' input.pathToTurbSim inputfile]); %add path to turbsim from input

if status ~= 0
    disp(cmdout)
    error(['Error in Turbsim execution file ' inputfile '. Aborting...'])   
else
    disp( ['.wnd ' inputfile ' is done!'])
end

%Delete the copied input files from the output folder
delete([input.pathToTurbSim Filename '.TimeSer'])
delete([input.pathToTurbSim Filename '_ConTurbSim' '.inp'])

%% Transform into .mat from .wnd and save it
[velocity, ~, ~, ~, ~, dz, dy, dt, ~, ~, SummVars] = readBLgrid([destination Filename '_ConTurbSim'  '.wnd']); % read the windfield .wnd
windfield = velocity2windfield(velocity,dz,dy,dt,SummVars); %#ok<*NASGU> % Convert to .mat
parsave('-v7.3',[input.TurbSimOut_dir Filename '_ConTurbSim' '.mat'],windfield)
disp([Filename '_ConTurbSim' ' wind field is converted to .mat'])


