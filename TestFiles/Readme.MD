# Test case 

The files provided here are meant to be used as a test to make sure everything works qwith the framework. In order to run the example:

1. Copy `InputParametersExample.m` to the *Main* folder 
2. Rename it to `InputParameters.m` and keep the original `InputParameters.m` file as a backup
3. Create a new folder inside the repository at the same level with *Main* folder called *OriginalWF*
4. Copy the *DTU10MW_Sh15_TI05_V08.mat* file in the newly created *OriginalWF* folder
3. Run `ViConDAR.m` file to test all modules.

Comments:
- To create a constrained wind field with Turbsim *Turbsim.exe* should be in the dedicated folder (look at `InputParametersExample.m` file)
- Python along with PyConturb package and its dependencies need to be installed in your computer in order to use PyConTurb
- Running PyConTurb through matlab (selected option in `InputParametersExample.m`)
