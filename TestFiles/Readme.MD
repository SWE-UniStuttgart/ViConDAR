# Test case 

The files provided here are meant to be used as a test to make sure everything works with the framework. In order to run the example:

1. Copy `InputParametersExample_fixed.m` or the `InputParametersExample_floating.m` to the *Main* folder. The former is an example of a forward scanning nacelle lidar application considering no dynamic motion of the lidar. The latter is an example of a forward scanning nacelle lidar application considering also the dynamic motions of the lidar, in this case created based on a floating wind turbine simulation using the openFAST tool from NREL. 
2. Rename it to `InputParameters.m` and keep the original `InputParameters.m` file as a backup
3. Create a new folder inside the repository at the same level with *Main* folder called *OriginalWF*
4. Copy the *TEST_Sh25_SD01_V15_TI05.mat* file for the fixed example and the *FLTEST_Sh20_SD01_V10_TI09.mat* for the floating example in the newly created *OriginalWF* folder
5. For the floating example create a new folder inside the repository at the same level with *Main* folder called *Floater_Dynamics* and copy the *FLTEST_Sh20_SD01_V10_TI09_dynamics.mat* into the newly created folder
6. Run `ViConDAR.m` file to test all modules.

Comments:
- To create a constrained wind field with Turbsim *Turbsim.exe* should be in the dedicated folder (look at `InputParametersExample.m` file)
- Be aware that the resampling of the lidar outputs with the resampling option may cause problems in the exectution of Turbsim in some cases. Refer to Turbsim input files for debugging. 
- Python along with PyConturb package and its dependencies need to be installed in your computer in order to use PyConTurb
- Running PyConTurb through matlab (selected option in `InputParametersExample.m`) requires python activated as a system variable (in Windows) and expects PyConTurb (along with its dependencies) to be installed in the base environment

If you don't have python in your system or PyConTurb installed in your base environment you can switch off flag *flag_obtain_Con_Pyconturb_Matlab* and switch on *flag_obtain_Con_Pyconturb_python*.
This will create a .csv input for the python wrapper `Wrapper_Read_CSV_Run_PycConTurb.py` located in *Functions* folder.
After this is created you need to run the wrapper script in python (e.g. through PyCharm or Spyder) in the PyConTurb environment. 
When the simulation is done, run ViConDAR again switching on only the flags *flag_obtain_Con_Pyconturb_python* and *flag_obtain_Con_Pyconturb_ConverToMat* to obtain the file in the global windfiled format.
Now you can run again ViConDAR using only the flags *flag_calculate_fullWF_statistics*, *flag_plot_lidar*, *flag_plot_WF_timeseries*, and *plot_fullWF_Slices* to finish testing all functionalities.

For the consideration of floating dynamics in the lidar simulation, the dynamics input time series must be provided in an additional dynamics inputs file (example: *FLTEST_Sh20_SD01_V10_TI09_dynamics.mat*)
All input channels (see example) must be available in the dynamics input file with same length. If dynamics channels are not available populate with zeros.
Here we need to add some description for the dynamic inputs, the structure and the requirements. 

You can check the lidar output files in the folder *LidarOutput*, the constrained wind fields in *ConstrainedWF\Turbsim* and *ConstrainedWF\PyConTurb* and the statistics in *ConstrainedWF\Statistics*.
For the batch post processing capabilities you can check the scripts in the folder *HelpfulStandAlone*.
