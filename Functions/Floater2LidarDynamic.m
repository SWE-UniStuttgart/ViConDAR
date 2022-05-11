%% Floater2LidarDynamic
% Function:Calculates translational displacement from 6DoF input
% considering offset between reference CS and Lidar CS, calculates
% translational velocity if calc_vel flag is set to true
%
%% Usage:
%
% [Output] = Floater2LidarDynamic(input,curFileInfo, Output)

%

%
%% Input:
%

%
%% Output:
%
%Output.x_I;
%Output.y_I;
%Output.z_I;
%
%
%% Modified:
%
%
%
%% ToDo:
%
%
%
%% Created:
%Moritz Gräfe 12.07.21

%% Code:
function [Output] = Floater2LidarDynamic(input, Output)


samplingfreq=Output.sim.dt;
%n=length(input.PatternY{1, 1}); %number of points in pattern
% LiDAR position in Floater CS
x_FL= input.x_FL;
y_FL= input.y_FL;
z_FL= input.z_FL;

Pitch= deg2rad(eval(input.channels.pitch))';
Roll=deg2rad(eval(input.channels.roll))';

Surge=eval(input.channels.x_trans_L)';
Sway= eval(input.channels.y_trans_L)';
Heave=eval(input.channels.z_trans_L)';

[x_trans,y_trans,z_trans] = arrayfun(@(Pitch, Roll) calc(Pitch, Roll, x_FL, y_FL, z_FL), Pitch, Roll);

    function [x_trans_Li,y_trans_Li,z_trans_Li]=calc(Pitch, Roll, x_FL, y_FL,z_FL)
        Yaw=0;
        
        % Yaw is a rotation around z-axis
        T_Yaw 	= [ cos(Yaw)   -sin(Yaw)	0;
            sin(Yaw) 	cos(Yaw)    0;
            0           0           1];
        
        % Pitch is a rotation around y-axis
        T_Pitch = [	cos(Pitch)	0           sin(Pitch);
            0         	1        	0;
            -sin(Pitch)	0        	cos(Pitch)];
        
        % Roll is a rotation around x-axis
        T_Roll  = [	1           0       	0;
            0       	cos(Roll)  -sin(Roll);
            0          	sin(Roll)	cos(Roll)];
   
        T       = T_Yaw*T_Pitch*T_Roll;

        x_trans_Li     = T(1,1)*x_FL + T(1,2)*y_FL + T(1,3)*z_FL;
        y_trans_Li     = T(2,1)*x_FL + T(2,2)*y_FL + T(2,3)*z_FL;
        z_trans_Li     = T(3,1)*x_FL + T(3,2)*y_FL + T(3,3)*z_FL;
   
    end

Output.x_trans_Li=((x_trans+ Surge)-x_FL);                            % - original pos % Lidar position in Inertial CS
Output.y_trans_Li=(y_trans+ Sway)-y_FL;                               % - original pos
Output.z_trans_Li=(z_trans+ Heave)-z_FL;    

Output.sim.TimeResults.x_trans_Li.Data={Output.x_trans_Li'};
Output.sim.TimeResults.y_trans_Li.Data={Output.y_trans_Li'};
Output.sim.TimeResults.z_trans_Li.Data={Output.z_trans_Li'};

% calculate tranlational velocities based on position or take velocity
% input time series

Output.x_vel=(diff(Output.x_trans_Li)/samplingfreq);
Output.y_vel=(diff(Output.y_trans_Li)/samplingfreq);
Output.z_vel=(diff(Output.z_trans_Li/samplingfreq));


Output.x_vel(:,end+1:length(Output.x_trans_Li))=nan;
Output.y_vel(:,end+1:length(Output.x_trans_Li))=nan;
Output.z_vel(:,end+1:length(Output.x_trans_Li))=nan;

Output.sim.TimeResults.x_vel.Data={Output.x_vel'};
Output.sim.TimeResults.y_vel.Data={Output.y_vel'};
Output.sim.TimeResults.z_vel.Data={Output.z_vel'};
    

end