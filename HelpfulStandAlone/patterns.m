%% Header
% Create different lidar patterns in cartesian coordinates. Origin is the
% WT hub center
% F. Costa
% University of Stuttgart, Stuttgart Wind Energy (SWE) 2019

%Example:
%     flag_square5          = 0;
%     flag_square9          = 1;
%     flag_circular         = 0;
%     center                = [0,0]; % center of the pattern (m)
%     
%     % Data of square patterns:
%     up                    = 30; % Z dimension of square pattern (m)
%     right                 = 30; % Y dimension of square pattern (m)
%     
%     
%     % Data of circular pattern:
%     radii                 = 20;    % radius of circle (m)
%     step_circular_pattern = 10;    % number of points along circle
%     
%     [Y,Z] = patterns (up,right,center,radii,step_circular_pattern,flag_square5, flag_square9, flag_circular);

%-------------------------------------------------------------------------------------------------------------------

function [Y,Z]=patterns (up,right, center,radii,step_circular_pattern,flag_square5, flag_square9, flag_circular)

        % Square pattern (5 points)
        if flag_square5==1
            
            Point1 = center;
            Point2 = [up,right];
            Point3 = [-up,right];
            Point4 = [up,-right];
            Point5 = [-up,-right];
            
            Y= [Point1(1),Point2(1),Point3(1),Point4(1),Point5(1)];
            Z= [Point1(2),Point2(2),Point3(2),Point4(2),Point5(2)];
        end
        
        % Square pattern (9 points)
        if flag_square9==1
            
            Point1 = center;
            Point2 = [up,right];
            Point3 = [-up,right];
            Point4 = [up,-right];
            Point5 = [-up,-right];
            Point6 = [0,-right];
            Point7 = [0,right];
            Point8 = [up,0];
            Point9 = [-up,0];
            
            Y = [Point1(1),Point2(1),Point3(1),Point4(1),Point5(1),Point6(1),Point7(1),Point8(1),Point9(1)];
            Z = [Point1(2),Point2(2),Point3(2),Point4(2),Point5(2),Point6(2),Point7(2),Point8(2),Point9(2)];
        end
        
        % Circular pattern:
        if flag_circular==1
            
            centers = center;
            radio   = radii;
                        
            viscircles(centers,radio);
            
            datos_circulo=findall(gca,'type','line');
            
            X2 = get(datos_circulo,'xdata');
            Y2 = get(datos_circulo,'ydata');
                        
            Y = X2{1}(1:step_circular_pattern:182);
            Z = Y2{1}(1:step_circular_pattern:182);
            hold on
            plot(Y,Z,'*b','LineWidth',2)            
        end       
end
