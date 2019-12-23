# Virtual Constrained turbulence and liDAR measurements (ViConDAR)

## Introduction

Modular framework allowing to scan a numerical wind field with a virtual lidar and create constrained turbulent wind fields for wind turbine simulations. 

The aim of the framework is to provide an open source base for such tests but also to be used as a platform where other modules (e.g. different velocity vector reconstruction methods) can be easily added and tested. More information on motivation and possible applications can be found in [1].

![flowchart](Pictures_repo/vicondar_overview.png)

## Code and folder structure

The framework is based around a wrapper script calling the different modules. The input/output definition is done through a function that gathers all the necessary parameters called *InputParameters*. In principal this should be the only interaction with the user. The different functionalities (e.g. performing virtual lidar measurements or constrianing the windfields) are activated by using the relevant flags defines in the inputs.  

### Main folder

#### InputParameters.m
The main IO interface of ViConDAR. Here you can switch on/off flags, define parameters and request plots

#### ViConDAR.m
Main wrapper script to be executed

### Functions Folder
All functions required for VIConDAR are gathered here. Any new additions should be added here.

### HelpfulStandAlone folder
This folder contains scripts used for batch processing and other small useful tools. These run as stand alone and do not interface directly with ViConDAR
### Other directories
The other directories are defined in the *InputParameters.m* and are created while the code is running. These directories include folders for lidar outputs, inputs to constrained turbulence generators as welll as outputs

## Run Test Case

The files provided in the TestCase folder are meant to be used as a test to make sure everything works qwith the framework. In order to run the example:

1. Copy *InputParametersExample.m* to the Main folder 
2. Rename it to *InputParameters.m* and keep the original *InputParameters.m* file as a backup
3. Create a new folder inside the repository at the same level with *Main* folder called *OriginalWF*
4. Copy the *DTU10MW_Sh15_TI05_V08.mat* file in the newly created *OriginalWF* folder
3. Run `ViConDAR.m` file to test all modules.

Comments:
- To create a constrained wind field with Turbsim Turbsim.exe should be in the dedicated folder (look at InputParametersExample.m file)
- Python along with PyConturb package and its dependencies need to be installed in your computer in order to use PyConTurb
- Running PyConTurb through matlab (selected option in InputParametersExample.m)
 
## Adding user modules


## Dependencies - Requirements
- Currently the constrained turbulence generators connected to ViConDAR are Turbsim (https://github.com/OpenFAST/openfast/tree/master/modules/turbsim) and PyConTurb (https://gitlab.windenergy.dtu.dk/rink/pyconturb). 
- Some functions regarding reading and writing .wnd files are created and distributed as open source by NREL.
- Python(>v3.5) has to be installed in the system. The path to the correct python.exe corresponding to the python environment including PyConTurb and its dependencies has to be explicitly stated in the InputParameters file.  
- Ther framework is developed and tested with Matlab 2015b 64bit. We have done some tests with Matlav 2018b 64bit but not all functionalities are guaranteed to work properly

## Resources
1. A numerical framework for constraining synthetic wind fields with lidar measurements for improved load simulations V.Pettas, F.Costa, J. Rinker, M. Kretschmer, P.W. Cheng
2. Rinker, J. M., “PyConTurb: an open-source constrained turbulence generator,” Journal of Physics: Conference Series, Vol.
1037, 2018, p. 062032. URL http://stacks.iop.org/1742-6596/1037/i=6/a=062032?key=crossref.e906e6eebbe4b3f8b1b8ad8e740718b1.
3. NWTC Information Portal (Alpha Versions).  https://nwtc.nrel.gov/Alphas. Last modified 14-June-2016 ;

## Citing
If you use ViConDAR for a publication, use the following citation:

A numerical framework for constraining synthetic wind fields with lidar measurements for improved load simulations V.Pettas, F.Costa, J. Rinker, M. Kretschmer, P.W. Cheng

