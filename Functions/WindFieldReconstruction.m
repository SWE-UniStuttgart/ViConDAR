function [VFinalTotal_Time_U,VFinalTotal_Time_V,VFinalTotal_Time_W] = WindFieldReconstruction(VFinalTotal_Time_LOS_vec, LOS_2_In_matrix, input, trajectory_forAng, ref_plane_dist,Y)
% Reconstructs the wind vector from the LiDAR  measurements based on the method
% specified by input.ReconstructionMethod
%
% Reconstruction options:
%    Simple: This reconstruction uses the v=0, w=0 assumption to reconstruct the u component of the wind speed for each pattern point
%    global_projection:  All measurement  from one Scan are used to reconstruct one estimation of u- and v-
%        component of the wind vector. This method requires at least 3
%        pattern point per scan. Based on https://elib.uni-stuttgart.de/handle/11682/8813
%
% TODO: add calculation of inflow angle
%
% V.Pettas/F.Costa/M.Gräfe
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2021

switch input.ReconstructionMethod
    case 'simple'
        for ind_LOS = 1:length(Y)
            for ind_slice = 1:size(VFinalTotal_Time_LOS_vec,2)
                VFinalTotal_Time_reconstr_vec = [ LOS_2_In_matrix{ind_LOS}(1,:); [0 0 0]; [0 0 0]] *...
                    [VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(1);VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(2) ;VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(3) ] ;
                
                VFinalTotal_Time_U{ind_LOS}(1,ind_slice) = VFinalTotal_Time_reconstr_vec(1); %#ok<*AGROW>
                VFinalTotal_Time_V{ind_LOS}(1,ind_slice) = VFinalTotal_Time_reconstr_vec(2);
                VFinalTotal_Time_W{ind_LOS}(1,ind_slice) = VFinalTotal_Time_reconstr_vec(3);
            end
        end
        
    case 'global_projection'
        for ind_slice = 1:size(VFinalTotal_Time_LOS_vec,2)
            for ind_LOS = 1:length(input.PatternY{1,1})
                VLOS(ind_LOS,1) = sqrt(VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(1)^2+VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(2)^2+VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice}(3)^2);
                A(ind_LOS,1 )   = ref_plane_dist/sqrt(trajectory_forAng(1, ind_LOS)^2+trajectory_forAng(1, ind_LOS)^2+ref_plane_dist^2);
                A(ind_LOS,2)    = trajectory_forAng(1,ind_LOS)/sqrt(trajectory_forAng(1, ind_LOS)^2+trajectory_forAng(1, ind_LOS)^2+ref_plane_dist^2);
            end
            S = pinv(A)*VLOS;
            for ind_LOS = 1:length(input.PatternY{1,1})
                VFinalTotal_Time_U{ind_LOS}(1,ind_slice) = S(1);
                VFinalTotal_Time_V{ind_LOS}(1,ind_slice) = S(2);
                VFinalTotal_Time_W{ind_LOS}(1,ind_slice) = 0;
            end
        end
        % Possible TODO: add calculation of inflow angle
end
end

