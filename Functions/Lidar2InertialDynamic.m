function [Dynamics] = Lidar2InertialDynamic(input, Dynamics)

% Converts the pattern points (x_L, y_L, z_L) from the lidar
% coordinate system to the earth fixed coordinte system (x_I, y_I, z_I)
% considering the roll, pitch and yaw rotational displacement
% V.Pettas/M.Gräfe
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021
%% Code
x_L = -Dynamics.x_L; % The minus is to account for the directeion of the coordinate system (x-axis pointing downwind)
y_L = Dynamics.y_L;
z_L = Dynamics.z_L;

channelnames = input.channels.channelnames;
Pitch = deg2rad(Dynamics.sim.channels_time.(channelnames{1})); % the input is in degrees
Roll  = deg2rad(Dynamics.sim.channels_time.(channelnames{2}));
Yaw   = deg2rad(Dynamics.sim.channels_time.(channelnames{3}));

[x_I,y_I,z_I] = arrayfun(@(Pitch, Roll,Yaw, x_L, y_L, z_L) calc(Pitch, Roll,Yaw, x_L, y_L, z_L), Pitch, Roll,Yaw, x_L, y_L, z_L);

% Remove offset from y_I, z_I channels. Why is this here in the bottom???
Dynamics.x_I = x_I;
Dynamics.y_I = y_I+input.Pos_LiDAR(1);
Dynamics.z_I = z_I+input.Pos_LiDAR(2);


    function [x_I,y_I,z_I] = calc(Pitch, Roll,Yaw, x_L, y_L,z_L)
        % Yaw displacement
        T_Yaw 	= [ cos(Yaw)   -sin(Yaw)	0;
            sin(Yaw) 	cos(Yaw)    0;
            0           0           1];
        % Pitch displacement
        T_Pitch = [	cos(Pitch)	0           sin(Pitch);
            0         	1        	0;
            -sin(Pitch)	0        	cos(Pitch)];
        % Roll displacement
        T_Roll  = [	1           0       	0;
            0       	cos(Roll)  -sin(Roll);
            0          	sin(Roll)	cos(Roll)];
        % Combination of rotational displacement
        T       = T_Yaw*T_Pitch*T_Roll;
        
        x_R     = T(1,1)*x_L + T(1,2)*y_L + T(1,3)*z_L;
        y_R     = T(2,1)*x_L + T(2,2)*y_L + T(2,3)*z_L;
        z_R     = T(3,1)*x_L + T(3,2)*y_L + T(3,3)*z_L;
        
        x_I     = -x_R;
        y_I     =  y_R;
        z_I     =  z_R;
    end


end
