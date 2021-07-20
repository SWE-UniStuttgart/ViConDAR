# Test case 

The files provided here are meant to be used as a test to make sure everything works qwith the framework. In order to run the example:

1. Copy `InputParametersExample.m` to the *Main* folder 
2. Rename it to `InputParameters.m` and keep the original `InputParameters.m` file as a backup
3. Create a new folder inside the repository at the same level with *Main* folder called *OriginalWF*
4. Copy the *TEST_Sh25_SD01_V15_TI05.mat* file in the newly created *OriginalWF* folder
3. Run `ViConDAR.m` file to test all modules.

Comments:
- To create a constrained wind field with Turbsim *Turbsim.exe* should be in the dedicated folder (look at `InputParametersExample.m` file)
- Python along with PyConturb package and its dependencies need to be installed in your computer in order to use PyConTurb
- Running PyConTurb through matlab (selected option in `InputParametersExample.m`) requires python activated as a system variable (in Windows) and expects PyConTurb (along with its dependencies) to be installed in the base environment

If you don't have python in your system or PyConTurb installed in your base environment you can switch off flag *flag_obtain_Con_Pyconturb_Matlab* and switch on *flag_obtain_Con_Pyconturb_python*.
This will create a .csv input for the python wrapper `Wrapper_Read_CSV_Run_PycConTurb.py` located in *Functions* folder.
After this is created you need to run the wrapper script in python (e.g. through PyCharm or Spyder) in the PyConTurb environment. 
When the simulation is done, run ViConDAR again switching on only the flags *flag_obtain_Con_Pyconturb_python* and *flag_obtain_Con_Pyconturb_ConverToMat* to obtain the file in the global windfiled format.
Now you can run again ViConDAR using only the flags *flag_calculate_fullWF_statistics*, *flag_plot_lidar*, *flag_plot_WF_timeseries*, and *plot_fullWF_Slices* to finish testing all functionalities.

You can check the lidar output files in the folder *LidarOutput*, the constrained wind fields in *ConstrainedWF\Turbsim* and *ConstrainedWF\PyConTurb* and the statistics in *ConstrainedWF\Statistics*.
For the batch post processing capabilities you can check the scripts in the folder *HelpfulStandAlone*.
