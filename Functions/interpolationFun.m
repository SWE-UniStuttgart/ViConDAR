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
                for iTSlice = 1:length(LOS_points.slicesAv) % For measured slices in the pattern
                    indLoopT =  (LOS_points.slices(i1,:)-1+LOS_points.slicesAv(iTSlice))*input.distanceSlices; % Distances where the measurements are focused  along the probe length
                    indLoopT2 = indLoopT;
                    indNEg = find(indLoopT<=0); % find negative, zeros or Nans
                    indNEg =  [indNEg find(isnan(indLoopT))]; % find negative or Nans
                    indNEg = [indNEg find(indLoopT>size(component,2)*input.distanceSlices)]; % findd points outside of the grid (we assume squared grid)
                    indLoopT2(indNEg) = [];
                    maxLoopInd = round(length(indLoopT)/2); % find the middle of the total slices
                    NansStart  = length(find(indNEg<maxLoopInd));
                    NansEnd    = length(indNEg)-NansStart;
                    indLoopT2=[nan(1,NansStart) indLoopT2 nan(1,NansEnd) ];                    
                    for ind_p=1:size(Y1,2)  % for points in the probe length
                            xq1 = Y1(1,ind_p);
                            xq2 = Z1(1,ind_p);
                            xq3 = indLoopT2;               
                            
                            % Interpolation
                            VFinalTotal_Time{i1}(ind_p,:)=interpn(gridy,input.slicesDistance,gridz,component,xq1,xq3,xq2); 
                    end
                end
       
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
        
        % Take the values at these points:
        if length(LOS_points.slicesAv)~=1
            VFinalTotal{i1} = squeeze(component(PoinInd{i1}(2,floor(length(PoinInd{i1})/2+1)),:,PoinInd{i1}(1,floor(length(PoinInd{i1})/2+1)))); % the floor is done to always take the middle point of consecuteive measurement ranges                     
            VFinalTotal{i1} = round(VFinalTotal{i1},5);
        else
            VFinalTotal{i1} = squeeze(component(PoinInd{i1}(2,1),:,PoinInd{i1}(1,1))); % the floor is done to always take the middle point of consecuteive measurement ranges
            VFinalTotal{i1} = round(VFinalTotal{i1},5);
        end
        

     % ####### Check the interpolation is working properly#####################               

    %                 busquedaY = find(ismember(gridy,Y1{dim1,dim2}(1,ind_p)) ); % look for coincidences in Y component
    %                 busquedaZ = find(ismember(gridz,Z1{dim1,dim2}(1,ind_p))); %#ok<*EFIND> % look for coincidences in  Z component
    %                 if isempty(busquedaY)%|| length(busquedaY) <(length(LOS_points.slicesAv))
    %                     DifY = gridy - Y1{dim1,dim2}(1,ind_p);
    %                     [~,ind_miny] = mink(abs(DifY),2,2);
    %                     ind_miny=sort(ind_miny);
    %                     Point_ind(1,:) = [(ind_miny(1)) (ind_miny(2))];
    %                     Point(1,:) = [gridy(ind_miny(1)) gridy(ind_miny(2))];
    %                 else
    %                     Point_ind(1,:) = [busquedaY,busquedaY];
    %                     Point(1,:)= [gridy(busquedaY) gridy(busquedaY)];
    %                 end
    %                 if isempty(busquedaZ)%|| %length(busquedaZ) <(length(LOS_points.slicesAv))
    %                     DifZ = gridz-Z1{dim1,dim2}(1,ind_p);
    %                     [~,ind_minz] = mink(abs(DifZ),2,2); % we get the two closest points (previous and following ones)
    %                     ind_minz=sort(ind_minz);
    %                     Point_ind(2,:) = [(ind_minz(1)) (ind_minz(2))];
    %                     Point(2,:) = [gridz(ind_minz(1)) gridz(ind_minz(2))];
    %                     
    %                 else
    %                     Point_ind(2,:) = (busquedaZ);
    %                     Point(2,:)= [gridy(busquedaZ) gridy(busquedaZ)];
    %                 end       
    %                 Point_ind(3,:)=input.focus_distances_index{ind_p};
    %                 Point(3,:) = input.focus_distances_new{ind_p};
    %                 point_to_interpolate{ind_p} = (combvec(Point_ind(1,:),Point_ind(3,:),Point_ind(2,:))');
    %                 point_to_interpolate_val{ind_p} = (combvec(Point(1,:),Point(3,:),Point(2,:))');
    %                 % Get velocity values
    %                 for ind_velo_vec=1:size(point_to_interpolate{ind_p},1)
    %                     point_to_interpolate_val{ind_p}(ind_velo_vec,4)= component(point_to_interpolate{ind_p}(ind_velo_vec,1),point_to_interpolate{ind_p}(ind_velo_vec,2),point_to_interpolate{ind_p}(ind_velo_vec,3));
    %                 end
    %                 % interpolate
    %                 point_to_interpolate_val{ind_p}=sortrows(point_to_interpolate_val{ind_p},2);
    %                 x1  = unique(point_to_interpolate_val{ind_p}(:,1)); % Y component
    %                 x2  = flipud(unique(point_to_interpolate_val{ind_p}(:,3))); % Z component
    %                 x3  = unique(point_to_interpolate_val{ind_p}(:,2)); % focus distance (X component)
    %                 Vel_val  = point_to_interpolate_val{ind_p}(:,4); % velocity vector

                    % cases for interpolation: In total there are 8 different
                    % cases to interpolate (have to implement the rest)

    %                 if isempty(busquedaY)==0  && isempty( busquedaZ)==0 && length(x3)~=1 % if Y and Z exist in the grid, but not x interpolate only in the x direction (along the focus distance)
    %                     Vel_mat{1,dim2}(1,ind_p) = interp1(x3, unique(Vel_val,'stable'),xq3,'linear');
    % %                 
    % %                 elseif isempty(busquedaY)==1  && isempty( busquedaZ)==0 
    % %                     Vel_mat{1,dim2}(1,ind_p) = interp2(x1,x3, unique(Vel_val,'stable'),xq1,xq3);
    % %                 
    % %                 elseif isempty( busquedaY)==0  && isempty( busquedaZ)==1
    % %                     Vel_mat{1,dim2}(1,ind_p) = interp2(x2,x3, unique(Vel_val,'stable'),xq2,xq3);
    % %                 
    % %                 elseif isempty(busquedaY)==1  && isempty(busquedaZ)==1 && length(x3)==1
    % %                     Vel_mat{1,dim2}(1,ind_p) = interp2(x1,x2, unique(Vel_val,'stable'),xq1,xq2);                     
    %                 
    %                 elseif isempty( busquedaY)==1  && isempty( busquedaZ)==1 && length(x3)~=1
    %                     V   = reshape(Vel_val,2,2,2);
    %                     Vel_mat{1,dim2}(1,ind_p)= interpn(x1,x2,x3,V,xq1,xq2,xq3,'linear');
    %                 end                                
    % ##########################################################################

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
            VFinalTotal{i1} = squeeze(component(PoinInd{i1}(2,1),:,PoinInd{i1}(1,1))); % the floor is done to always take the middle point of consecuteive measurement ranges
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
        if length(LOS_points.slicesAv) ~= 1
            VFinalTotal_Time{i1} = mean(VFinalTotal_TimeInt2,'omitnan');% Change it for a gaussian mean!!! Averaging columns which contain all the volume averaging ppins in the LOS
        else
            VFinalTotal_Time{i1} = VFinalTotal_TimeInt2;
        end
        clear VFinalTotal_TimeInt
    end
end
