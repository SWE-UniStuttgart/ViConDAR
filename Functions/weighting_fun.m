function VFinalTotal_Time=weighting_fun(input,LOS_points,VFinalTotal_TimeInt2)
if length(LOS_points.slicesAv) ~= 1
    if strcmpi(input.flag_probe_weighting,"mean")
        for ind_mean=1:size(VFinalTotal_TimeInt2,2)
            VFinalTotal_Time{ind_mean} = mean(VFinalTotal_TimeInt2{ind_mean},'omitnan');% Change it for a gaussian mean!!! Averaging columns which contain all the volume averaging ppins in the LOS
        end
    
    elseif strcmpi(input.flag_probe_weighting,"gaussian")
        % Introducing Gaussina weights in the performance of probe volume:
        % First we create the weights:
        for i=1:size(VFinalTotal_TimeInt2,2)
            VFinalTotal_TimeInt3=VFinalTotal_TimeInt2(:,i);
            
%             VFinalTotal_TimeInt_noNAN=VFinalTotal_TimeInt3(~isnan(VFinalTotal_TimeInt3)); % remove nans
%             weights = linspace (-length(VFinalTotal_TimeInt_noNAN),length(VFinalTotal_TimeInt_noNAN),size(VFinalTotal_TimeInt_noNAN,1));
            
            % For a given fwhm (fwhm = 2*Rayleigh length) we calculate the normal distribution:
            fwhm           = 2*input.distance_av_space;
            focus_distance = input.ref_plane_dist;
            offset = 100;
            distan = linspace(focus_distance-offset,focus_distance+offset,1e4);
            
            %Sigma and fwhm are simply related as follows:
            sigma          = fwhm/2.355;
            gaussian       = (1/(sigma*sqrt(2*pi)))*exp(-0.5*((distan-focus_distance)/sigma).^2);
            
            %Find the points above the fwhm
            half_max       = (min(gaussian) + max(gaussian)) / 2;
            probe_distance = distan(gaussian >= half_max); % distances above the fwhm
            
            % Distances at slices queried by the lidar may not match the distances wihin the probe volume. We take the closest points. Other method might be interpolate        
            gaussian_factor = zeros(1, length(input.focus_distances));% preallocation
            % Round to the millimeter
            distan         = round(distan,4);
            probe_distance = round(probe_distance,4);
            for dist=1:length(input.focus_distances) % for all the queried points within the probe volume                
                query_point = find (distan==input.focus_distances(dist),1);
                if isempty (query_point) % if the queried point is not in the distances vector created to perform the gaussian weights, take the closest
                    [~,ind] =min(abs(probe_distance-input.focus_distances(dist)));
                    gaussian_factor(1,dist) = gaussian (distan==probe_distance(ind)); % gaussian factors for the distances queried by the lidar within the probe length                                       
                else
                    gaussian_factor(1,dist) = gaussian (query_point);                    
                end
            end
            
            % probability_weights = gaussian_factor*(distan(2)-distan(1));
            % performing weighted mean
            VFinalTotal_TimeInt3_NoNans                  = isnan(VFinalTotal_TimeInt3); %finding nans
            gaussian_factor(VFinalTotal_TimeInt3_NoNans) = nan;
            VFinalTotal_Time(:,i)                        = sum(gaussian_factor'.*VFinalTotal_TimeInt3,'omitnan')/sum(gaussian_factor,'omitnan');
                       
            % Using cumulative distribution to perform probabilities
%             cumulative_probability = cumsum(gaussian)*(distan(2)-distan(1));
%             for ind_cum=1:length(input.focus_distances)
%                 [~,pos_cum]=min(abs(distan-input.focus_distances(ind_cum)));
%                 index_cum{ind_cum}=[pos_cum,pos_cum+1];
%             end
%             
%             for ind_cum1=1:length(index_cum)
%                 cum_prob(1,ind_cum1) =  (distan(index_cum{ind_cum1}(2))*cumulative_probability(index_cum{ind_cum1}(2)))-(distan(index_cum{ind_cum1}(1))*cumulative_probability(index_cum{ind_cum1}(1)));
%             end
% 
%             % performing weighted mean with cumulative distribution:
%             VFinalTotal_Time2(:,i) = sum(cum_prob'.*VFinalTotal_TimeInt3,'omitnan');

        end
    end
else
    VFinalTotal_Time = VFinalTotal_TimeInt2;
end


