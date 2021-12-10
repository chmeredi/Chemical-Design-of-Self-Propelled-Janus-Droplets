classdef profile_class6
   % Flowprofile class
   % Input: 4 matrixes x,y,u,v containing x and y coordinates and x and y
   % directions of flow on those coordinates
   % To prepare data for use of this object class:
   % 'Load("directory\nameofdata")' % Keep in mind data should have format
   % x,y,u,v where all are m x n matrices of the same size
   % 
   %
   % Methods:
   %    - initialization method
   %    - rotate_profile: rotates by angle radians. Keeps track of angle
   %    - get_angle: outputs rotation angle with respect to original frame
   %    - obj_size: outputs size of grid
   %    - translate_profile: translates by transX and transY pixels. Keeps track of translation
   %    - get_location: outputs location with respect to original frame
   %    - plot profile: plots the flow profile
   %    - reform_grid : plots profile on a different grid new_x, new_y
   %    - scale_profile: scales flow vectors in profile by a factor scale
   %    - substract_profile: if two profiles have the same grid: substracts the new profile from object
   %    - rotate_back: rotates profile back to original orientation
   %    - translate_back: translates profile back to original orientation
   %    - obj=reform_grid(obj,new_x,new_y)
   %    - obj=get_new_grid(listofprofileclasses)
   %    - obj=average_profiles(listofprofileclasses)
   %    - obj=average_profiles_samegrid(listofprofileclasses)
   %    - obj=average_profiles_fast(listofprofileclasses)
   %    - plot_poles(obj)
   %    - obj_angle=get_pole_orient_cartesian(obj)
   %    - obj_polar=cart2pol_poles(obj)
   %    - draw_angle_polar(obj)
   %    - output=get_border_distance(obj,centerx,centery)
   %    - obj=crop_profile(obj,centerx,centery,radius)
   %
   %
   % Created: on 08-08-16 by Pepijn Moerman
   % Last modified: 11-10-16 by Pepijn Moerman
   properties
       x        % nxn matrix containing x coordinates of gridpoints
       y        % nxn matrix containing y coordinates of gridpoints
       u        % nxn matrix containing x coordinates of flow at gridpoints x,y
       v        % nxn matrix containing y coordinates of flow at gridpoints x,y
       ang      % if profile has been rotated: provides final rotation angle with original profile
       X        % if profile has been tranlated: provides final x-translation
       Y        % if profile has been tranlated: provides final y-translation
   end
   properties
       rotation_Num = 0;        % number of times the profile was rotated
       translation_Num = 0;     % number of times the profile was translated
       all_angles = [];         % 1xN matrix containing all subsequent rotation angles
       all_positions = [];      % 2xN matrix containing all subsequent translations
   end
   properties (Dependent)
       obj_spacing       % spacing between gridpoints in flow profile
       obj_size          % gridsize of type [# x gridpoints, # y gridpoints]
       obj_angle         % in case of quadrupolar flow profile: angle of profile. Calculated in cartesian coordinates)
       obj_polar         % nx2 matrix containing all velocities in polar coordintes [theta,rho]
   end
   methods 
       function obj=profile_class6(x,y,u,v)
           obj.x = x;
           obj.y = y;
           obj.u = u;
           obj.v = v;
           
       end
       function obj=rotate_profile(obj,angle)
            %Rotates flow velocities and grid by a an angle 'angle'
            
            obj_size=get_object_size(obj);
            %Creates empy matrices for rotated flows
            unew=NaN(obj_size);
            vnew=NaN(obj_size);
            xnew=NaN(obj_size);
            ynew=NaN(obj_size);
            
            %Rotates the flow velocities
            for i=1:obj_size(2)
                for j=1:obj_size(1)
                    ui=obj.u(j,i);
                    vi=obj.v(j,i);
                    %create coordinate vector
                    coord=[ui,vi];
                    %create rotation matrix and multiply with coordinate vector
                    b = [cos(angle)  -sin(angle) ; sin(angle)  cos(angle)] * coord';
                    unew(j,i)=b(1);
                    vnew(j,i)=b(2);
                end
            end
    
            %Rotate coordinate system flow profile
            for i=1:obj_size(2)
                for j=1:obj_size(1)
                    xi=obj.x(j,i);
                    yi=obj.y(j,i);
                    %create coordinate vector
                    coord=[xi,yi];
                    %create rotation matrix and multiply with coordinate vector
                    b = [cos(angle)  -sin(angle) ; sin(angle)  cos(angle)] * coord';
                    xnew(j,i)=b(1);
                    ynew(j,i)=b(2);
                end
            end
            %Save rotation angle
            obj.rotation_Num=obj.rotation_Num+1;
            obj.all_angles(obj.rotation_Num)=angle;
            obj.ang=sum(obj.all_angles);
            obj.x=xnew;
            obj.y=ynew;
            obj.u=unew;
            obj.v=vnew;
       end
       function rotatedby=get_angle(obj)
           %outputs the angle by which the object has been rotated
           rotatedby=obj.ang;
       end
       function obj_size=get_object_size(obj)
           %Gets the size of the grid of object
           obj_size=[length(obj.x(:,1)),length(obj.x(1,:))];
       end
       function obj_spacing= get_object_spacing(obj)
           x_bound=length(obj.x(1,:));
           y_bound=length(obj.y(:,1));
           obj_spacing=NaN;
           %Gets the spacing of the grid
           for i=1:x_bound-1
               for j=1:y_bound
                   obj_spacing=sqrt((obj.x(j,i+1)-obj.x(j,i))^2+(obj.y(j,i+1)-obj.y(j,i))^2);
                   if isnan(obj_spacing)==0;
                       break
                   end
               end
               if isnan(obj_spacing)==0;
                   break
               end
           end
           
       end
       function location=get_location(obj)
           %returns X and Y coordinate of translation
           location = [obj.X,obj.Y];
           
       end
       function obj=translate_profile(obj,transX,transY)
           %Translates the profile by over a vector transX,transY
           obj.translation_Num = obj.translation_Num + 1;               % Counts number of translation
           obj.x=obj.x+transX;                                          % Translates x coordinates
           obj.y=obj.y+transY;                                          % Translates y coordinates
           obj.all_positions(1:2,obj.translation_Num)=[transX,transY];  % appends list of all translations
           obj.X=sum(obj.all_positions(1,:));                           % saves current x position with respect to original
           obj.Y=sum(obj.all_positions(2,:));                           % saves current y position with respect to original
           
       end
       function plot_profile(obj,boundaries)
           %Plots the profile
           %boundaries must be a 1x4 matrix that indicates the plotrange
           %[xmin xmax ymin ymax]
           
           figure
           hold on
           box on
           if ~nargin<2
               axis(boundaries)
           end
           quiver(obj.x,obj.y,obj.u,obj.v,2,'b');
           
           hold off
       end
       function obj=scale_profile(obj,scale)
           %Scales the flow velocities of the profile by a factor 'scale'
           obj.u=obj.u*scale;
           obj.v=obj.v*scale;
       end
       function obj=translate_back(obj)
           %Translates the object back to its original position
           location=get_location(obj);
           obj=translate_profile(obj,-location(1),-location(2));
       end
       function obj=rotate_back(obj)
           %Rotates the object back to its original position
           angle=get_angle(obj);
           obj=rotate_profile(obj,-angle);
       end
       function obj=substract_profile(obj,other_profile)
           %substracts other_profile off of obj
           %other profile must be of type profile_class
           %still to fill in.
           %grids of obj and other_profile must be the same
           
           %assert that grids have the right shape
           if obj.x~=other_profile.x | obj.y~=other_profile.y
               error('Grids of two profiles are not the same. Consider using reform_grid to make sure both grids have the same shape')
           end
           
           %substract flow profiles
           obj.u=obj.u-other_profile.u;
           obj.v=obj.v-other_profile.v;
           
       end
       function obj=reform_grid(obj,new_x,new_y)
           % Reshapes the grid on which a flow profile is plotted to a new
           % grid of arbitrary shape 'new_x,new_y'
           % Assumes both grids have isotropic spacing
           
           % Get spacing of old and new grid and take maximum spacing
           obj_size=get_object_size(obj);
           obj_spacing_old=get_object_spacing(obj);
           obj_spacing_new=abs(new_x(1,2)-new_x(1,1));
           obj_spacing=max(obj_spacing_old,obj_spacing_new);
           
           % Initialize empty matricees
           new_u=NaN(size(new_x));
           new_v=NaN(size(new_x));
           %note: obj_size is [height,width] so [y,x]
           for i=1:length(new_x(1,:))  % for all x coordinates on the new grid
               for j=1:length(new_x(:,1)) % for all y coordinates on the new grid
                   count=1; %Initialize counter and empty matrices
                   distance=[];
                   Uadd=[];
                   Vadd=[];
                   weight=[];
                   for g=1:obj_size(2) %for all x coordinates in old frame
                       for h=1:obj_size(1) % for all y coordinates in the old frame
                           %calculate distance between point on old grid
                           %g,h and new grid i,j
                           distance(h,g)=sqrt((new_x(j,i)-obj.x(h,g))^2+(new_y(j,i)-obj.y(h,g))^2);
                           if distance(h,g) <= obj_spacing
                               weight(count)=obj_spacing-distance(h,g); %Save distance to weigh later
                               % Calculate flow contribution of each point weighed by distance
                               Uadd(count)=obj.u(h,g);
                               Vadd(count)=obj.v(h,g);
                               % Raise the counter by one, to prevent overwriting
                               count=count+1;
                           end
                       end
                   end
                   % Fill new grid with averaged flow vectors by summing weighted contributions
                   weighted_u=NaN(count-1,1);
                   weighted_v=NaN(count-1,1);
                   max_weight=sum(weight(:));
                   for q=1:count-1
                       weighted_u(q)=Uadd(q)*weight(q)/max_weight;
                       weighted_v(q)=Vadd(q)*weight(q)/max_weight;
                   end
                   new_u(j,i)=sum(weighted_u(:));
                   new_v(j,i)=sum(weighted_v(:));
                   
               end
           end
       % Write new grid as output
       obj.x=new_x;
       obj.y=new_y;
       obj.u=new_u;
       obj.v=new_v;
       end
       function obj=get_new_grid(listofprofileclasses)
           % Here listofprofileclasses is a single object of type profile class 
           % or a matrix containing n objects of type profile_class
           % Returns an empty grid large enough to fit all grids in the list
           % Grid has largest spacing present in list
           
           % Obtain dimensions and location of input grids
           size_list=size(listofprofileclasses);
           max_x_value=NaN(size_list(2),1);
           min_x_value=NaN(size_list(2),1);
           max_y_value=NaN(size_list(2),1);
           min_y_value=NaN(size_list(2),1);
           spacing=NaN(length(listofprofileclasses),1);
           for n=1:length(listofprofileclasses)
               profile_n=listofprofileclasses(n);
               spacing(n)=get_object_spacing(profile_n);
               max_x_value(n)=max(max(profile_n.x));
               min_x_value(n)=min(min(profile_n.x));
               max_y_value(n)=max(max(profile_n.x));
               min_y_value(n)=min(min(profile_n.x));
           end
           
           % Fetch largest spacing and extreme x and y points
           new_spacing=max(spacing);
           new_spacing(isnan(new_spacing(:)))=[];
           x_lowest=min(min_x_value);
           x_lowest(isnan(x_lowest(:)))=[];
           x_highest=max(max_x_value);
           x_highest(isnan(x_highest(:)))=[];
           y_lowest=min(min_y_value);
           y_lowest(isnan(y_lowest(:)))=[];
           y_highest=max(max_y_value);
           y_highest(isnan(y_highest(:)))=[];
           
           % Calculate number of steps to cover all space with new_spacing
           nsteps_x=round((x_highest(1)-x_lowest(1))/new_spacing(1));
           nsteps_y=round((y_highest(1)-y_lowest(1))/new_spacing(1));
           
           % Create the new grid
           new_x=NaN(nsteps_y,nsteps_x);
           new_y=NaN(nsteps_y,nsteps_x);
           new_u=NaN(nsteps_y,nsteps_x);
           new_v=NaN(nsteps_y,nsteps_x);
           
           % Fill the x coordinates of the grid
           for i=1:nsteps_x
               new_x(:,i)=x_lowest+new_spacing*(i-1);
           end
           
           for j=1:nsteps_y
               new_y(j,:)=y_lowest+new_spacing*(j-1);
           end
           
           obj=profile_class6(new_x,new_y,new_u,new_v);
           
       end
       function obj=average_profiles(listofprofileclasses)
           % Here listofprofileclasses is an array containing n objects of
           % type profile_class
           
           'This may take a while'
           
           % Gets the new grid
           new_grid=get_new_grid(listofprofileclasses);
           gridsize=size(new_grid.x);
           new_u=NaN(gridsize);
           new_v=NaN(gridsize);
           gridspace=get_object_spacing(new_grid);
           
           %For all points on the grid: that is, for all x coordinates...
           for i=1:gridsize(2)
               %...and for all y coordinates
               for j=1:gridsize(1)
                   %Uadd and Vadd will contain all velocities u,v from neighbouring points g,h in frame k to be averaged
                   %over in the point i,j
                   Uadd=[];
                   Vadd=[];
                   distance=[];
                   weight=[];
                   %and start a counter
                   count=1;
                   %at each point on the grid: loop through all frames
                   for n=1:length(listofprofileclasses)
                       %and in each frame, go through all points on the frame...
                       %...in x
                       profile_n=listofprofileclasses(n);
                       size_profile_n=size(profile_n.x);
                       for h=1:size_profile_n(2)
                           %...and in y
                           for g=1:size_profile_n(1)
                               %calculate the distance between this point g,h in this
                               %frame k and the point on the grid i,j
                               distance(h,g)=sqrt((new_grid.x(j,i)-profile_n.x(h,g))^2+(new_grid.y(j,i)-profile_n.y(h,g))^2);
                               %if the distance is smaller than the gridspacing
                               if distance(h,g) < gridspace
                                   weight(count)=gridspace-distance(h,g); %Save distance to weigh later
                                   %then we want to remember it:
                                   %so we add the value of u for this point to the list Uadd
                                   %and the value of v to the list Vadd
                                   Uadd(count)=profile_n.u(h,g);
                                   Vadd(count)=profile_n.v(h,g);
                                   %raise the counter by one, to prevent overwriting
                                   count=count+1;
                               end
                           end
                       end
                   end
                   weighted_u=NaN(count-1,1);
                   weighted_v=NaN(count-1,1);
                   max_weight=sum(weight(:));
                   for q=1:count-1
                       weighted_u(q,1)=Uadd(q)*weight(q)/max_weight;
                       weighted_v(q,1)=Vadd(q)*weight(q)/max_weight;
                   end
                   weighted_u(isnan(weighted_u(:,1)))=[];
                   weighted_v(isnan(weighted_v(:,1)))=[];
                   new_u(j,i)=sum(weighted_u(:));
                   new_v(j,i)=sum(weighted_v(:));
               end
           end
           
           obj=profile_class(new_grid.x,new_grid.y,new_u,new_v);
           
       end 
       function obj=average_profiles_samegrid(listofprofileclasses)
           % Here listofprofileclasses is an array containing n objects of
           % type profile_class
           % Assumes all grids in listofprofileclasses lie are the same
           
           x_new=listofprofileclasses(1).x;
           y_new=listofprofileclasses(1).y;
           grid_size=size(x_new);
           u_new=NaN(grid_size);
           v_new=NaN(grid_size);
           
           for i=1:grid_size(2)
               for j=1:grid_size(1)
                   Uadd=NaN(length(listofprofileclasses),1);
                   Vadd=NaN(length(listofprofileclasses),1);
                   for n=1:length(listofprofileclasses)
                       profile_n=listofprofileclasses(n);
                       Uadd(n,1)=profile_n.u(j,i);
                       Uadd(isnan(Uadd(:,1)))=[];
                       Vadd(n,1)=profile_n.v(j,i);      
                       Vadd(isnan(Vadd(:,1)))=[];
                   end
                   u_new(j,i)=mean(Uadd(:,1));
                   v_new(j,i)=mean(Vadd(:,1));
               end
           end
           
           obj=profile_class6(x_new,y_new,u_new,v_new);
       end
       function obj=average_profiles_fast(listofprofileclasses)
           % Here listofprofileclasses is an array containing n objects of
           % type profile_class
           % Returns the average of all itmes in listofprofileclasses
           
           new_grid=get_new_grid(listofprofileclasses);
           new_list=[];
           gridsize=size(new_grid.x);
           if gridsize(1)>200 ^ gridsize(2)>200
               disp('Warning: at least 1 dimension of gridsize is exceptionally large. Be aware that this will cause the program to take a long time to finish.')
           end
           for n=1:length(listofprofileclasses)
               if ~mod(n/10,1)
                disp([num2str(length(listofprofileclasses)-n) ' frames to go.']);
               end
               item=listofprofileclasses(n);
               item=reform_grid(item,new_grid.x,new_grid.y);
               new_list=[new_list item];
           end
           obj=average_profiles_samegrid(new_list);
       end
       function plot_poles(obj)
           % Make 1D matrix of all values of u and v
           obj_size=get_object_size(obj);
           u_stretch=NaN(obj_size(1)*obj_size(2),1);
           v_stretch=NaN(obj_size(1)*obj_size(2),1);
           for i=1:obj_size(2)
               u_stretch(1+(i-1)*obj_size(1):i*obj_size(1),1)=obj.u(:,i);
               v_stretch(1+(i-1)*obj_size(1):i*obj_size(1),1)=obj.v(:,i);
           end
          figure;
          plot(u_stretch,v_stretch,'b.')
       end
       function obj_angle=get_pole_orient_cartesian(obj)
           % Fitting procedure to obtain the orientation of the flow
           % profile
           % Takes a cone with angle_width and looks for points in that
           % cone
           % Rotates by angle_step and repeat
           % Returns angle with maximum points in cone
           
           % Get initial search values
           angle_width=pi/4;
           angle_step=pi/20;
           counter=2;
           round_angle=NaN(7,1);
           best_angle=NaN;
           for i=1:7
               round_angle(i,1)=1000*i;
           end
           
           % Get polar angles
           obj_polar=cart2pol_poles(obj);
           
           while abs(round_angle(counter)-round_angle(counter-1))>0.1 && counter<7
               % Initialize counters
               count_old=0;
               nsteps=round((pi/2)/angle_step);
               
               % For all values in the anglerange pi/2
               for a=0:nsteps
                   angle=a*angle_step;
                   count=0;
                   % Find all points in the four cones related to that value
                   for n=1:length(obj_polar(:,1))
                       % If the point n is in the primary cone
                       if obj_polar(n,1)>(angle-angle_width/2) && obj_polar(n,1)<(angle+angle_width/2)
                           count=count+1;
                           % If the point n is in the opposite cone
                           if obj_polar(n,1)>(pi+angle-angle_width/2) && obj_polar(n,1)<(pi+angle+angle_width/2)
                               count=count+1;
                           end
                           % If the point n is in the left side cone
                           if obj_polar(n,1)>(pi/2+angle-angle_width/2) && obj_polar(n,1)<(pi/2+angle+angle_width/2)
                               count=count+1;
                           end
                           % If the point n is in the right side cone
                           if obj_polar(n,1)>(-pi/2+angle-angle_width/2) && obj_polar(n,1)<(-pi/2+angle+angle_width/2)
                               count=count+1;
                           end
                       end
                   end
                   % Keep the best value for count and angle
                   if count>count_old
                       best_angle=a*angle_step;
                       count_old=count;
                   end
               end
               % Remember the best angle of this iteration
               % Add up the counter
               % and make the searching range smaller
               round_angle(counter-1)=best_angle;
               counter=counter+1;
               angle_width=angle_width/2;
               angle_step=angle_step/2;
           end
           
           % Return the best angle value
           obj_angle=round_angle(counter-2);
       end          
       function obj_polar=cart2pol_poles(obj)
           % Calculates polar coordinate version of u and v
           
           % Make 1D matrix of all values of u and v
           obj_size=get_object_size(obj);
           u_stretch=NaN(obj_size(1)*obj_size(2),1);
           v_stretch=NaN(obj_size(1)*obj_size(2),1);
           for i=1:obj_size(2)
               u_stretch(1+(i-1)*obj_size(1):i*obj_size(1),1)=obj.u(:,i);
               v_stretch(1+(i-1)*obj_size(1):i*obj_size(1),1)=obj.v(:,i);
           end
           
           [theta,rho]=cart2pol(u_stretch(:,1),v_stretch(:,1));
           for i=1:length(theta(:,1))
               if theta(i,1)<0
                   theta(i,1)=2*pi-abs(theta(i,1));
               end
           end
           obj_polar(:,1:2)=[theta rho];
       end
       function draw_angle_polar(obj)
           % Plots u versus v and draws two diagonals separating the 4
           % highest point density areas in the plot
           
           % Make 1D matrix of all values of u and v
           obj_size=get_object_size(obj);
           u_stretch=NaN(obj_size(1)*obj_size(2),1);
           v_stretch=NaN(obj_size(1)*obj_size(2),1);
           for i=1:obj_size(2)
               u_stretch(1+(i-1)*obj_size(1):i*obj_size(1),1)=obj.u(:,i);
               v_stretch(1+(i-1)*obj_size(1):i*obj_size(1),1)=obj.v(:,i);
           end
           
           % Get extrema of polar_plot
           min_value=1.2*min(u_stretch);
           max_value=1.2*max(u_stretch);
           min_y=1.2*min(v_stretch);
           max_y=1.2*max(v_stretch);
           
           % Get angle
           obj_angle=get_pole_orient_polar(obj);
           
           % Calculate lines
           line1=@(x) tan(obj_angle+pi/4)*x;
           line2=@(x) tan(obj_angle-pi/4)*x;
           
           % Draw the plot and the lines
           figure;
           hold on
           fplot(line1,[min_value max_value],'k-')
           fplot(line2,[min_value max_value],'k-')
           plot(u_stretch,v_stretch,'b.')
           set(gca,'Ylim',[min_y,max_y],'Xlim', [min_value,max_value]);
           hold off
           
           
       end
       function output=get_border_distance(obj,centerx,centery)
           % Finds the edge of the profile from a central point
           % centerx,centery and outputs minimal distance from that
           % point to the center
           
           % Find edge
           distance(1,1)=abs(centerx-max(max(obj.x)));
           distance(2,1)=abs(centerx-min(min(obj.x)));
           distance(3,1)=abs(centery-max(max(obj.y)));
           distance(4,1)=abs(centery-min(min(obj.y)));
           min_dist=min(distance);
           
           % Output the minimal distance from point to edge
           output=min_dist;
           
       end
       function obj=crop_profile(obj,centerx,centery,radius)
           % Takes a flow profile obj of N by N points, a central point and
           % a radius.
           % Creates a circle of radius around central point
           % centerx,centery
           % Changes flow vectors outside the circle to NaN
           % Outputs the new flow profile obj.
           
           % Get gridsize
           object_size=get_object_size(obj);
           
           % Change all points outside circle to NaN
           for i=1:object_size(1)
               for j=1:object_size(2)
                   dist=sqrt((centerx-obj.x(i,j))^2+(centery-obj.y(i,j))^2);
                   if dist>radius
                       obj.x(i,j)=NaN;
                       obj.y(i,j)=NaN;
                       obj.u(i,j)=NaN;
                       obj.v(i,j)=NaN;
                   end
               end
           end
                 
       end
       function obj_angle=get_pole_orient_quadrupole(obj)
           % Returns the angle of a quadrupolar flow profile
           % If flow profile is not quadrupolar, returns NaN
           
           % Get polar coordinates of flow
           obj_polar=cart2pol_poles(obj);
           max_radius=max(obj_polar(:,2));
           cutoff=max_radius*0.4;
           
           % Take all angles of flows with a length larger than cutoff
           data=[];
           counter=0;
           
           for i=1:length(obj_polar(:,1))
               if obj_polar(i,2)>cutoff
                   counter=counter+1;
                   data(counter,1)=obj_polar(i,1);
               end
           end
           
           
           
           % Try orientations and calculate points in the four poles
           
           % Initialize counters and search bounderies
           angle_step=pi/200;
           count_old=0;
           best_angle=NaN;
           boundaries_all=NaN(4,2);
           boundaries_all(:,1)=pi/8;
           boundaries_all(:,2)=-pi/8;
           boundaries_all(2,:)=boundaries_all(2,:)+pi/2;
           boundaries_all(3,:)=boundaries_all(3,:)+pi;
           boundaries_all(4,:)=boundaries_all(4,:)+3*pi/2;
           
           for alpha=0:angle_step:pi/2;
               % Reset counter
               negativecheck=zeros(4,1);
               count1=0;
               count2=0;
               count3=0;
               count4=0;
               boundaries(:,:)=boundaries_all(:,:)+alpha;
               % Remove negatives from boundaries matrix
               for bb=1:4
                   if boundaries(bb,1)<0
                       boundaries(bb,1)=2*pi+boundaries(bb,1);
                       negativecheck(bb,1)=1;
                   end
                   if boundaries(bb,2)<0
                       boundaries(bb,2)=2*pi+boundaries(bb,2);
                       negativecheck(bb,1)=1;
                   end
               end
               % For all points in data
               for n=1:length(data(:,1))
                   
                   % If point is in regime 1
                   if negativecheck(1,1)==1
                       if data(n,1)<boundaries(1,1) || data(n,1)>boundaries(1,2);
                           count1=count1+1;
                       end
                   else
                       if data(n,1)<boundaries(1,1) && data(n,1)>boundaries(1,2)
                           count1=count1+1;
                       end
                   end
                   % If point is in regime 2
                   if data(n,1)<boundaries(2,1) && data(n,1)>boundaries(2,2)
                       count2=count2+1;
                   end
                   % If point is in regime 3
                   if data(n,1)<boundaries(3,1) && data(n,1)>boundaries(3,2)
                       count3=count3+1;
                   end
                   % If point is in regime 4
                   if negativecheck(4,1)==1
                       if data(n,1)<boundaries(4,1) || data(n,1)>boundaries(4,2);
                           count4=count4+1;
                       end
                   else
                       if data(n,1)<boundaries(4,1) && data(n,1)>boundaries(4,2)
                           count4=count4+1;
                       end
                   end
               end
               count_tot=count1+count2+count3+count4;
               % Keep the best value for count and angle
               if count_tot>count_old
                   best_angle=alpha;
                   count_old=count_tot;
                   count1_best=count1;
                   count2_best=count2;
                   count3_best=count3;
                   count4_best=count4;
               end
               
           end
          
           if count_old<0.6*length(data(:,1)) || count1_best< 0.1*0.4*length(data(:,1)) || count2_best< 0.1*0.4*length(data(:,1)) || count3_best< 0.1*0.4*length(data(:,1)) || count4_best< 0.1*0.4*length(data(:,1));
               best_angle=NaN;
               'It looks like this frame is not a quadrupole'
               
           end
           
           % Rotate profile by calculated angle
           new_profile=rotate_profile(obj,-best_angle);
           
           % Get all elemants in the positive y-half of the profile
           all_pos_y=new_profile.y>0;
           all_v_there=new_profile.v(all_pos_y);
           all_v_there(isnan(all_v_there(:,1)),:)=[];
           mean_speed_there=mean(all_v_there);
           if mean_speed_there<0
               obj_angle=best_angle-pi/2;
           else
               obj_angle=best_angle;
           end
           

       end
       function decay=get_radial_decay(obj,centerx,centery)
           % For each gridpoints, calculates the radial distance to center
           % and the magnitude of the velocity vector
           % Plots the velocity versus the distance
           obj_size=get_object_size(obj);
           counter=0;
           all_distances=NaN(obj_size(1)*obj_size(2),1);
           all_speeds=NaN(obj_size(1)*obj_size(2),1);
           for i=1:obj_size(1)
               for j=1:obj_size(2)
                   counter=counter+1;
                   speed=sqrt(obj.u(j,i)^2+obj.v(j,i)^2);
                   distance=sqrt((obj.x(j,i)-centerx)^2+(obj.y(j,i)-centery)^2);
                   all_distances(counter,1)=distance;
                   all_speeds(counter,1)=speed;
               end
           end
           
           decay=[all_distances(:,1) all_speeds(:,1)];           
       end
       function decay=get_radial_decay_polar(obj,centerx,centery)
           % For each gridpoints, calculates the radial distance to center
           % and the magnitude of the radial component of the velocity vector
           % Plots the velocity versus the distance
           
           
           obj_size=get_object_size(obj);
           counter=0;
           %Create empty matrices with distance and speeds
           all_distances=NaN(obj_size(1)*obj_size(2),1);
           all_v_r=NaN(obj_size(1)*obj_size(2),1);
           all_v_theta=NaN(obj_size(1)*obj_size(2),1);
           %Loop through all points in x,and y position
           for i=1:obj_size(2)
               for j=1:obj_size(1)
                   %count length of list
                   counter=counter+1;
                   %calculate polar components of speed
                   speed=sqrt(obj.u(j,i)^2+obj.v(j,i)^2);
                   angle_s=atan2(obj.v(j,i),obj.u(j,i));
                   %calculate absolute distance
                   distance=sqrt((obj.x(j,i)-centerx)^2+(obj.y(j,i)-centery)^2);
                   %calculate angle of position with drop
                   angle=atan2(obj.y(j,i)-centery,obj.x(j,i)-centerx);
                   %radial component of speed is cos(theta)*totalspeed
                   %with theta is angle of speed - angle of position
                   %(angle-angle_s)
                   v_r=cos(angle-angle_s)*speed;
                   v_theta=sin(angle-angle_s)*speed;
                   all_distances(counter,1)=distance;
                   all_v_r(counter,1)=v_r;
                   all_v_theta(counter,1)=v_theta;
               end
           end
           
           decay=[all_distances(:,1) all_v_r(:,1) all_v_theta(:,1)];           
       end
       function file=plot_radial_decay(obj,centerx,centery)
           % Plots the radial decay of a flow profile, binned with
           % errorbars
           
           % Calculate [radius, flow magnitude] matrix
           decay=get_radial_decay(obj,centerx,centery);
           
           % Assert the list is not too small to do binning
           if length(decay(:,1))<100
               'Object contains very few gridpoints. Modify binsize by hand to make this work. Or consider not binning at all'
           end
           
           npoints=100;
           min_rad=min(decay(:,1));
           max_rad=max(decay(:,1));
           binsize=(max_rad-min_rad)/npoints;
           binned_decay=NaN(100,1);
           binned_radius=NaN(100,1);
           
           for i=1:npoints
               bin_bound_low=min_rad+(i-1)*binsize;
               bin_bound_up=min_rad+i*binsize;
               magnitude_to_bin=[];
               radius_to_bin=[];
               for j=1:length(decay(:,1))
                   if decay(j,1)>bin_bound_low && decay(j,1)<bin_bound_up
                       magnitude_to_bin=[magnitude_to_bin decay(j,2)];
                       radius_to_bin=[radius_to_bin decay(j,1)];
                   end
               end
               
               if length(magnitude_to_bin)>0
                   binned_decay(i,1)=mean(magnitude_to_bin(1,:));
                   binned_decay(i,2)=std(magnitude_to_bin(1,:));
                   binned_radius(i,1)=mean(radius_to_bin(1,:));
                   binned_radius(i,2)=std(radius_to_bin(1,:));
               end
           end
           
           figure
           errorbar(binned_radius(:,1),binned_decay(:,1),binned_decay(:,2),'ko')
           
           % Returns a matrix containing, binned radius with error bars and
           % binned flow speed with error bars
           file=[binned_radius,binned_decay];
           
       end
       function obj=crop_profile_to_size(obj,centerx,centery,cropsize)
           % Crops profile to circular shape of size distance
           % The center is centerx, centery
           % All profile points outside of crop circle are turned to NaN
           
           % Get gridsize
           object_size=get_object_size(obj);
           
           % Find edge
           distance(1,1)=abs(centerx-max(max(obj.x)));
           distance(2,1)=abs(centerx-min(min(obj.x)));
           distance(3,1)=abs(centery-max(max(obj.y)));
           distance(4,1)=abs(centery-min(min(obj.y)));
           boundary=min(distance);
           
           % Get size
           radius=min([cropsize boundary]);
           
           % Change all points outside circle to NaN
           for i=1:object_size(1)
               for j=1:object_size(2)
                   dist=sqrt((centerx-obj.x(i,j))^2+(centery-obj.y(i,j))^2);
                   if dist>radius
                       obj.x(i,j)=NaN;
                       obj.y(i,j)=NaN;
                       obj.u(i,j)=NaN;
                       obj.v(i,j)=NaN;
                   end
               end
           end
       end
       function obj_angle=get_pole_orient_dipole(obj)
           % Returns the angle of a dipolar flow profile
           % If flow profile is not dipolar, returns NaN
           
           % Get polar coordinates of flow
           obj_polar=cart2pol_poles(obj);
           max_radius=max(obj_polar(:,2));
           cutoff=max_radius*0.4;
           
           % Take all angles of flows with a length larger than cutoff
           data=[];
           counter=0;
           
           for i=1:length(obj_polar(:,1))
               if obj_polar(i,2)>cutoff
                   counter=counter+1;
                   data(counter,1)=obj_polar(i,1);
               end
           end
           
           
           if length(data)>0
               % Try orientations and calculate points in the four poles
               
               % Initialize counters and search bounderies
               angle_step=pi/200;
               count_old=0;
               best_angle=NaN;
               boundaries_all=NaN(1,2);
               boundaries_all(:,1)=pi/16;
               boundaries_all(:,2)=-pi/16;
               
               for alpha=0:angle_step:2*pi;
                   % Reset counter
                   negativecheck=0;
                   count1=0;
                   boundaries(:,:)=boundaries_all(:,:)+alpha;
                   % Remove negatives from boundaries matrix
                   if boundaries(1,1)<0
                       boundaries(1,1)=2*pi+boundaries(1,1);
                       negativecheck=1;
                   end
                   if boundaries(1,2)<0
                       boundaries(1,2)=2*pi+boundaries(1,2);
                       negativecheck=1;
                   end
                   
                   % For all points in data
                   for n=1:length(data(:,1))
                       
                       % If point is in regime 1
                       if negativecheck==1
                           if data(n,1)<boundaries(1,1) || data(n,1)>boundaries(1,2);
                               count1=count1+1;
                           end
                       else
                           if data(n,1)<boundaries(1,1) && data(n,1)>boundaries(1,2)
                               count1=count1+1;
                           end
                       end
                   end
                   
                   % Keep the best value for count and angle
                   if count1>count_old
                       best_angle=alpha;
                       count_old=count1;
                   end
                   
               end
               
               if count_old<0.3*length(data(:,1))
                   best_angle=NaN;
                   'It looks like this frame has little dipolar charachter.'
                   
               end
               
               obj_angle=best_angle;
           else
               'Flow profile cannot be obtained. Probably droplet is only partly in view'
               obj_angle=NaN;
           end
           
       end
       function vector=get_vector_product(obj,centerx,centery)
           % For flow profile obj with droplet at position centerx, centery
           % Calculates for each point on the flow, the out product of
           % distance and flow velocity vector. Then sums the out products
           % Negative vector product means rotation one way, positive menas
           % rotation the other way (possibly)
           
           % Get object size
           obj_size=get_object_size(obj);
           prod=NaN(obj_size);
           
           % Loop through all x coordinates
           for i=1:obj_size(2)
               % Loop through all y coordinates
               for j=1:obj_size(1)
                   % Get coordinates and velocities
                   xi=obj.x(j,i)-centerx;
                   yi=obj.y(j,i)-centery;
                   ui=obj.u(j,i);
                   vi=obj.v(j,i);
                   a=sqrt(xi^2+yi^2);
                   b=sqrt(ui^2+vi^2);
                   phia=atan2(xi,yi);
                   phib=atan2(ui,vi);
                   theta=phia-phib;
                   if theta<-pi
                       theta=theta+2*pi;
                   elseif theta > pi;
                       theta=theta-2*pi;
                   end
                                      
                   % Calculate out product
                   prod(j,i)=b*sin(theta)/a^2; %For real out product 1=a
                   
               end
           end
           
           % Remove NaN from matrix
           prod(isnan(prod))=[];
           
           % Sum numbers to obtain total out product
           vector=sum(sum(prod));
           
                     
       end
   end
end
