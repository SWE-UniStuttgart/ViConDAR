function [Dynamics] = LidarTranslations(input, Dynamics)

% Extracts translational displacement and velocity channels from
% dynamics input and combines with lidar parameters to get the position, velocity,
% and displacement from the original position of the lidar system.
%
% V.Pettas/M.Gräfe
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021
%
%% Code:

channelnames=input.channels.channelnames;
% translational displacement Lidar
Dynamics.x_trans_Li = Dynamics.sim.channels_time.(channelnames{7});
Dynamics.y_trans_Li = Dynamics.sim.channels_time.(channelnames{8});
Dynamics.z_trans_Li = Dynamics.sim.channels_time.(channelnames{9});

% translational velocity of Lidar
Dynamics.x_vel = Dynamics.sim.channels_time.(channelnames{4});
Dynamics.y_vel = Dynamics.sim.channels_time.(channelnames{5});
Dynamics.z_vel = Dynamics.sim.channels_time.(channelnames{6});

%Lidar Position in Inertial CS
Dynamics.x_posL_I = Dynamics.x_trans_Li;
Dynamics.y_posL_I = Dynamics.y_trans_Li+input.Pos_LiDAR(1);
Dynamics.z_posL_I = Dynamics.z_trans_Li+input.Pos_LiDAR(2);

% Add translation of Lidar to x_I, y_I, z_I
Dynamics.x_I = Dynamics.x_I - Dynamics.x_trans_Li; % The minus is to account for the directeion of the coordinate system (x-axis pointing downwind)
Dynamics.y_I = Dynamics.y_I + Dynamics.y_trans_Li;
Dynamics.z_I = Dynamics.z_I + Dynamics.z_trans_Li;

end