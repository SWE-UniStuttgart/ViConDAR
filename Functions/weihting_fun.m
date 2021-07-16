function [VFinalTotal_Time]=weihting_fun(input,LOS_points,VFinalTotal_TimeInt2)
if length(LOS_points.slicesAv) ~= 1
    if strcmpi(input.flag_probe_weighting,"mean")
     VFinalTotal_Time{i1} = mean(VFinalTotal_TimeInt2,'omitnan');% Change it for a gaussian mean!!! Averaging columns which contain all the volume averaging ppins in the LOS
    % Introducing Gaussina weights in the performance of probe volume:
    % First we create the weights:
    elseif strcmpi(input.flag_probe_weighting,"gaussian")
    for i=1:size(VFinalTotal_TimeInt2,2)
        VFinalTotal_TimeInt3=VFinalTotal_TimeInt2(:,i);
        VFinalTotal_TimeInt_noNAN=VFinalTotal_TimeInt3(~isnan(VFinalTotal_TimeInt3)); % remove nans
        weights = linspace (-length(VFinalTotal_TimeInt_noNAN),length(VFinalTotal_TimeInt_noNAN),size(VFinalTotal_TimeInt_noNAN,1));
        distance
        pdf = fitdist(VFinalTotal_TimeInt_noNAN,'Normal'); %fiting to a normal distribution
        
        weights_gauss = normpdf(weights,0,pdf.sigma);
        
        % performing weighted mean
        VFinalTotal_Time{i1}(:,i) = sum(weights_gauss'.*VFinalTotal_TimeInt_noNAN)/sum(weights_gauss');
        
    end
    end
else
    VFinalTotal_Time{i1} = VFinalTotal_TimeInt2;
end