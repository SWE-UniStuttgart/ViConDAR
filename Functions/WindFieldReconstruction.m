function [VFinalTotal_Time_U,VFinalTotal_Time_V,VFinalTotal_Time_W] = WindFieldReconstruction(VFinalTotal_Time_LOS_vec, LOS_2_In_matrix, input,Y)
% Reconstructs the wind vector from the LiDAR  measurements based on the method
% specified by input.ReconstructionMethod
%
% Reconstruction options:
%    Simple: This reconstruction uses the v=0, w=0 assumption to reconstruct the u component of the wind speed for each pattern point
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
end
end

