%% Header
%
% Running pyconturb from matlab calling a python script translating .csv files
% to inputs for pyconturb and running it. Alternatively this can be done only
% in python by getting the name list and creting in loop in python 
%
% V.Pettas/F.Costa 
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019
 
function RunPyconturb_Con_WF(input,filename)

% For this step is is needed to have python >3.7 installed and added the path awith pyconturb package added.
% "python <fullpathtomyscript.py> <sys.argv1>... <sys.argv5>" arg1: contc file
% arg2: gridy file arg3: gridz file arg4: variables arg5: full path savename
[status,cmdout] = system([input.PythonExe_path ' ../Functions/Read_CSV_Run_PycConTurb.py ' input.PyconturbInput_dir 'con_tc_' filename '.csv ' input.PyconturbInput_dir 'GridY_' filename '.csv ' input.PyconturbInput_dir 'GridZ_' filename '.csv ' input.PyconturbInput_dir 'Variables_' filename '.csv ' input.PyconturbOut_dir filename '_ConPyconturb.csv']);
if status ~= 0
    disp(cmdout)
    error(['Error in Pyconturb execution file ' inputfile '. Aborting...'])   
end
% Transform from .csv into .mat
Pyconturb_CSV2MAT(input,filename)
disp([filename ' PyConTurb constrained windfield is done.'])
