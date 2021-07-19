function VFinalTotal_Time=weighting_fun(input,LOS_points,VFinalTotal_TimeInt2)
if length(LOS_points.slicesAv) ~= 1
    if strcmpi(input.flag_probe_weighting,"mean")
     VFinalTotal_Time = mean(VFinalTotal_TimeInt2,'omitnan');% Change it for a gaussian mean!!! Averaging columns which contain all the volume averaging ppins in the LOS
    
    
    elseif strcmpi(input.flag_probe_weighting,"gaussian")
        % Introducing Gaussina weights in the performance of probe volume:
        % First we create the weights:
        for i=1:size(VFinalTotal_TimeInt2,2)
            VFinalTotal_TimeInt3=VFinalTotal_TimeInt2(:,i);
%             f=find(isnan(VFinalTotal_TimeInt3)); %finding nans
%             VFinalTotal_TimeInt_noNAN=VFinalTotal_TimeInt3(~isnan(VFinalTotal_TimeInt3)); % remove nans
%             weights = linspace (-length(VFinalTotal_TimeInt_noNAN),length(VFinalTotal_TimeInt_noNAN),size(VFinalTotal_TimeInt_noNAN,1));
            
            % For a given fwhm ( fwhm = 2*Rayleigh length) we calculate the normal distribution:
            fwhm   = 2*input.distance_av_space;
            focus_distance   = input.ref_plane_dist;
            distan = linspace(-500+focus_distance-input.distance_av_space,500+focus_distance+input.distance_av_space,1e4);
            sigma=fwhm/2.355;
            gaussian = (1/(sigma*sqrt(2*pi)))*exp(-0.5*((distan-focus_distance)/sigma).^2);
            half_max=(min(gaussian) + max(gaussian)) / 2;

            index1 = find(gaussian >= half_max);
            probe_distance= distan(index1); % distances above the fwhm
            
            
            % Distances at slices queried by the lidar may not match the distances wihin the probe volume. We take the closest points. Other method might be interpolate
            
            for dist=1:length(input.focus_distances) % for all the queried points within the probe volume
                query_point = find (distan==input.focus_distances(dist),1);
                if isempty (query_point) % if the queried point is not in the distances vector created to perform the gaussian weights, take the closest
                    [~,ind] =min(abs(probe_distance-input.focus_distances(dist)));
                    gaussian_factor(1,dist) = gaussian (distan==probe_distance(ind));
                else
                    gaussian_factor(1,dist) = gaussian (distan==probe_distance(querypoint));
                end
                
            end


            
%             gaussian(f)=nan;
            probability_weights = gaussian_factor*(distan(2)-distan(1));
            
            
            
            
            

            % performing weighted mean
            VFinalTotal_Time(:,i) = sum(gaussian_factor'.*VFinalTotal_TimeInt3,'omitnan')/sum(gaussian_factor,'omitnan');
            VFinalTotal_Time2(:,i) = sum(probability_weights'.*VFinalTotal_TimeInt3,'omitnan');

        end
    end
else
    VFinalTotal_Time = VFinalTotal_TimeInt2;
end