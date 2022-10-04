function [VFinalTotal_Time_LOS_vec] = Weighting_Probe(VFinalTotal_Time_LOS_vec, Y)

% Adding noise to the measurements to imitate non-modelled noise in real measurements
% Currently gaussian white noise is added to assuming equal noise levels to all frequencies 
%   add ectr inputs needed for gaussian and peak detection 
%   add flag for changing weighting and peak detection methods

% add case falg differentiation
for ind_LOS = 1:length(Y)
    % the outpit has uvw components  in 1,2,3 dimension
    for ind_slice = 1:size(VFinalTotal_Time_LOS_vec,2)
        VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice} = mean(VFinalTotal_Time_LOS_vec{ind_LOS,ind_slice},2,'omitnan');
    end
end

end

