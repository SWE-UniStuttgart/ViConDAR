%% Header
%
% Interpolate to get all points (even when we are not in the grid) taking
% into account LOS and position of LiDAR.
%
% V.Pettas/F.Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019


function [VFinalTotal,VFinalTotal_Time,Y1,Z1] = interpolationFun(input,component,LOS_points,gridy,gridz,fullTime,windfield,dt,type_interpolation_2)

if input.interpolation_slices==1  %obsolete it shouldn't be used... Fix it later???? Currently we use the closest grid point and the nearest time slice
    for i1 = 1:size(LOS_points.slices,1) % loop over the points of pattern
        Y1 = LOS_points.Coor{i1}(1,:);
        Z1 = LOS_points.Coor{i1}(2,:);
        for ind_slice=1:size(LOS_points.slices,2)
            for iTSlice = 1:length(LOS_points.slicesAv) % For measured slices in the pattern

                indLoopT =  (LOS_points.slices(i1,ind_slice)+LOS_points.slicesAv(iTSlice))*input.distanceSlices  ;% [m] Distances where the measurements are focused  along the probe length
                indLoopT2 = indLoopT;
                indNEg = find(indLoopT<=0); % find negative, zeros or Nans
                indNEg =  [indNEg find(isnan(indLoopT))]; % find negative or Nans
                indNEg = [indNEg find(indLoopT>size(component,2)*input.distanceSlices)]; % findd points outside of the grid (we assume squared grid)
                indLoopT2(indNEg) = [];
                maxLoopInd = round(length(indLoopT)/2); % find the middle of the total slices
                NansStart  = length(find(indNEg<maxLoopInd));
                NansEnd    = length(indNEg)-NansStart;
                indLoopT2=[nan(1,NansStart) indLoopT2 nan(1,NansEnd) ];                    

                % Query points for interpolation:
                xq1 = LOS_points.Coor{i1}(1,iTSlice);
                xq2 = LOS_points.Coor{i1}(2,iTSlice);
                xq3 = indLoopT2;               
                % Interpolation
                VFinalTotal_TimeInt2{i1}(iTSlice,ind_slice)=interpn(gridz,input.slicesDistance,gridy,component,xq2,xq3,xq1); 
            end
        end
    end    
    % Apply weighting function:
    VFinalTotal_Time = weighting_fun(input,LOS_points,VFinalTotal_TimeInt2);
    
    % For the complete WF, insead of taking the closest, interpolate along time (x axis):
    for ind_s=1:size(LOS_points.Coor,2)
        Y2=LOS_points.Coor{ind_s}(1,:);
        Z2=LOS_points.Coor{ind_s}(2,:);
        Mid_Y2=Y2(floor(length(Y2)/2)+1); % find the point in the middle of the vector
        Mid_Z2=Z2(floor(length(Z2)/2)+1);
        VFinalTotal{ind_s} = interpn(gridz,input.slicesDistance,gridy,component,Mid_Z2,input.slicesDistance,Mid_Y2);
    end


else %if you don't interpolate get the closest point
    for i1 = 1:size(LOS_points.slices,1) % loop over the points of pattern
        Y1 = LOS_points.Coor{i1}(1,:);
        Z1 = LOS_points.Coor{i1}(2,:);
        Y1 = round(Y1,5); % Round because if not, the interpolation gives NaN's or incorrect results
        Z1 = round(Z1,5);
        
        busquedaY = find(ismember(gridy,Y1) ); % look for coincidences in Y component
        busquedaZ = find(ismember(gridz,Z1)); %#ok<*EFIND> % look for coincidences in  Z component
        if isempty (busquedaY) || length(busquedaY) <(length(LOS_points.slicesAv)) %check if all measured points are grid points
            for i = 1:length(Y1)    % find the closest points in the grid and use these
                DifY = gridy-Y1(i);
                [~,indF] = min(abs(DifY));
                Point_ind{i1}(1,i) = gridy(indF);
                PoinInd{i1}(1,i) = indF;
            end
        else
            Point_ind{i1}(1,:) = gridy(busquedaY);
            PoinInd{i1}(1,:) = busquedaY;
        end
        
        if isempty (busquedaZ) || length(busquedaZ) <(length(LOS_points.slicesAv)) %check if all measured points are grid points
            for i = 1:length(Y1)
                DifZ = gridz-Z1(i);
                [~,indF] = min(abs(DifZ));
                Point_ind{i1}(2,i) = gridz(indF);
                PoinInd{i1}(2,i) = indF;
            end
        else
            Point_ind{i1}(2,:) = gridz(busquedaZ); % points in the grid in meters matching our points (nearest, not exactly the value of the trajectory!!!)
            PoinInd{i1}(2,:) = busquedaZ; %indices of point matching in the grid
        end
        % Now take the values at this point
        if length(LOS_points.slicesAv)~=1
            VFinalTotal{i1} = squeeze(component(PoinInd{i1}(2,floor(length(LOS_points.slicesAv)/2+1)),:,PoinInd{i1}(1,floor(length(LOS_points.slicesAv)/2+1)))); % the floor is done to always take the middle point of consecuteive measurement ranges
            VFinalTotal{i1} = round(VFinalTotal{i1},5);
        else
            VFinalTotal{i1} = squeeze(component(PoinInd{i1}(2,1),:,PoinInd{i1}(1,1))); 
            VFinalTotal{i1} = round(VFinalTotal{i1},5);
        end
        
        for iTSlice = 1:length(LOS_points.slicesAv)% loop over the ranges of each point of the pattern and average them
            indLoopT =  (LOS_points.slices(i1,:)+LOS_points.slicesAv(iTSlice));
            indLoopT2 = indLoopT;
            indNEg = find(indLoopT<=0); % find negative, zeros or Nans
            indNEg =  [indNEg find(isnan(indLoopT))]; % find negative or Nans
            indNEg = [indNEg find(indLoopT>size(component,2))]; % findd points outside of the grid
            indLoopT2(indNEg) = [];
            VFinalTotal_TimeInt{iTSlice} = squeeze(component(PoinInd{i1}(2,iTSlice), indLoopT2 ,PoinInd{i1}(1,iTSlice)));    %
            
            % find how many nan points to add in the beginning or in the
            % end to keep length consistent
            maxLoopInd = round(length(indLoopT)/2); % find the middle of the total slices
            NansStart  = length(find(indNEg<maxLoopInd));
            NansEnd    = length(indNEg)-NansStart;
            VFinalTotal_TimeInt2(iTSlice,:) = [nan(1,NansStart) VFinalTotal_TimeInt{iTSlice} nan(1,NansEnd) ];
            
        end
        
        % Apply weighting function:
        VFinalTotal_TimeInt3{i1} =VFinalTotal_TimeInt2;
        VFinalTotal_Time = weighting_fun(input,LOS_points,VFinalTotal_TimeInt3);
             
%         if length(LOS_points.slicesAv) ~= 1
%             VFinalTotal_Time{i1} = mean(VFinalTotal_TimeInt2,'omitnan');% Change it for a gaussian mean!!! Averaging columns which contain all the volume averaging ppins in the LOS
%         else
%             VFinalTotal_Time{i1} = VFinalTotal_TimeInt2;
%         end

        clear VFinalTotal_TimeInt
    end
end
