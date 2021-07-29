%% Header
%
% Interpolate to get all points (even when we are not in the grid) taking
% into account LOS and position of LiDAR. 
%
% V.Pettas/F.Costa 
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019


function [VFinalTotal,VFinalTotal_Time,Y1,Z1] = interpolationFun(input,component,LOS_points,gridy,gridz,fullTime,dt,type_interpolation_2)

if input.interpolation_slices==1  %obsolete it shouldn't be used... Fix it later???? Currently we use the closest grid point and the nearest time slice
    for dim2=1:size (LOS_points.Coor,2)  
        for dim1=1:size (LOS_points.Coor,1)
            Y1{dim1,dim2} = LOS_points.Coor{dim1,dim2}(1,:);
            Z1{dim1,dim2} = LOS_points.Coor{dim1,dim2}(2,:);
            for ind_p=1:size(Y1{dim1,dim2},2)
                busquedaY = find(ismember(gridy,Y1{dim1,dim2}(1,ind_p)) ); % look for coincidences in Y component
                busquedaZ = find(ismember(gridz,Z1{dim1,dim2}(1,ind_p))); %#ok<*EFIND> % look for coincidences in  Z component
                if isempty(busquedaY)|| length(busquedaY) <(length(LOS_points.slicesAv))
                    DifY = gridy-Y1{dim1,dim2}(1,ind_p); 
                    [~,ind_miny] = mink(abs(DifY),2,2);
                    ind_miny=sort(ind_miny);
                    PointFm{dim1,dim2}(1,:) = [gridy(ind_miny(1)) gridy(ind_miny(2))];
                    
                else
                    PointFm{dim1,dim2}(1,:) = gridy(busquedaY);
                end
                if isempty(busquedaZ)|| length(busquedaZ) <(length(LOS_points.slicesAv))
                    DifZ = gridz-Z1{dim1,dim2}(1,ind_p); 
                    [~,ind_minz] = mink(abs(DifZ),2,2); % we get the two closest points (previous and following ones)
                    ind_minz=sort(ind_minz);
                    PointFm{dim1,dim2}(2,:) = [gridz(ind_minz(1)) gridz(ind_minz(2))];
                else
                    PointFm{dim1,dim2}(2,:) = gridz(busquedaZ);
                end
                
                % Combine points to interpolate them.
                n=1;
%                 for indfocus=1:size(input.focus_distances_new,2)
                    for ix=1:length(input.focus_distances_new{1})
                        for iy = 1:size(PointFm{dim1,dim2}(1,:),2) 
                            for iz=1:size(PointFm{dim1,dim2}(1,:),2) 
                                Point_to_interp{ind_p,n}=[vertcat(input.focus_distances_new{ind_p}(1,ix),PointFm{1,dim2}(1,iy),PointFm{1,dim2}(2,iz))];
                                n=n+1;
                            end
                        end
                    
                    end
            end
        end
    end
    
    
    
    
    
    %     for slice = 1:1:fullTime/dt+1
%         for i1 = 1:size(LOS_points.slices,1)
            
%             % Values with LOS of LiDAR:
%             Y1 = LOS_points.Coor{i1}(1,:);
%             Z1 = LOS_points.Coor{i1}(2,:);
%             Y1 = round(Y1,5); % Round because if not, the interpolation gives NaN's or incorrect results
%             Z1 = round(Z1,5);
%             if max(Y1)>max(gridy) || max(Z1)>max(gridz) || min(Y1)<min(gridy) || min(Z1)<min(gridz) % check if point are inside the grid
%                 % Do nothing
%             else
%                 if slice<=(fullTime/dt) % Delete the last slice to avoid errors in the process
%                     %% Space Interpolation.
%                     component1=squeeze(component(:,slice,:)); % Remove dimension=1                    
%                     VFinalTotal{i1}(:,slice)=interp2(gridy,gridz,component1,Y1,Z1,type_interpolation_2); %#ok<*AGROW>
%                     VFinalTotal{i1}(:,slice)=round(VFinalTotal{i1}(:,slice),2);
%                 end
%             end
%         end
%     end
    
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
                PointFm{i1}(1,i) = gridy(indF);
                PoinInd{i1}(1,i) = indF;
            end
        else
            PointFm{i1}(1,:) = gridy(busquedaY);
            PoinInd{i1}(1,:) = busquedaY;
        end
        
        if isempty (busquedaZ) || length(busquedaZ) <(length(LOS_points.slicesAv)) %check if all measured points are grid points
            for i = 1:length(Y1)
                DifZ = gridz-Z1(i);
                [~,indF] = min(abs(DifZ));
                PointFm{i1}(2,i) = gridz(indF);
                PoinInd{i1}(2,i) = indF;
            end
        else
            PointFm{i1}(2,:) = gridz(busquedaZ); % points in the grid in meters matching our points (nearest, not exactly the value of the trajectory!!!)
            PoinInd{i1}(2,:) = busquedaZ; %indices of point matching in the grid
        end
        % Now take the values at this point
        if length(LOS_points.slicesAv)~=1
            VFinalTotal{i1} = squeeze(component(PoinInd{i1}(2,floor(length(LOS_points.slicesAv)/2+1)),:,PoinInd{i1}(1,floor(length(LOS_points.slicesAv)/2+1)))); % the floor is done to always take the middle point of consecuteive measurement ranges
            VFinalTotal{i1} = round(VFinalTotal{i1},5);
        else
            VFinalTotal{i1} = squeeze(component(PoinInd{i1}(2,1),:,PoinInd{i1}(1,1))); % the floor is done to always take the middle point of consecuteive measurement ranges
            VFinalTotal{i1} = round(VFinalTotal{i1},5);
        end
        
        for iTSlice = 1:length(LOS_points.slicesAv)% loop over the ranges of each point of the pattern and average them
            indLoopT =  LOS_points.slices(i1,:)+LOS_points.slicesAv(iTSlice);
            
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
        if length(LOS_points.slicesAv) ~= 1            
            VFinalTotal_Time{i1} = mean(VFinalTotal_TimeInt2,'omitnan');% Change it for a gaussian mean!!! Averaging columns which contain all the volume averaging ppins in the LOS
        else
            VFinalTotal_Time{i1} = VFinalTotal_TimeInt2;
        end
        clear VFinalTotal_TimeInt
    end
end
