classdef PIVdroplets < droplets
    %PIVdroplets class is a sublcass of droplets
    %Incorporates functionalities of profile_class object and droplet
    %object
    %Allows further analysis of flow profiles obtained through PIVlab
    %Extra properties
    %   - Flows
    %   - average_flow
    %   - average_flow_rot
    %   - average_flow_trans
    %Extra methods
    %   - obj=add_PIVdata(obj)
    %   - obj=average_profile(obj,selection)
    %   - obj=average_profile_trans(obj)
    %   - obj=average_profile_rot(obj)
    %   - obj=substract_rot_trans(obj)
    %   - plot_flow(obj)
    %   - plot_flow_trans(obj)
    %   - plot_flow_rot(obj)
    %   - obj=substract_background_flow(obj)
    %   - plot_trace_on_movie_flow(obj,selection);
    %   - plot_flow_vs_distance_rot(obj)
    %
    %Uses external codes
    %   - profile_class6 and all contained within
    %Created on 11-01-18 by Pepijn Moerman
    %Last modified: 25-01-18 by Pepijn Moerman
    %               25-03-20 by Pepijn Moerman
    % - input of newer version of PIVlab now accepted
    % - produces contour correctly for newer version of pivlab
    % - ONLY works when in pivlab frames are split as 1-2, 3-4, NOT 1-2, 2-3!
    % - can deal with flow profiles around 2 particles as well as only 1
    %
    %How to use
    %Take the following steps to produce average flow profile:
    %       1. Create in ImageJ an 8 bit .tif stack with only the droplet of
    %       interest visibile entirely in every frame
    %       2. Create a folder with a list of individual .tif files with
    %       as the tiff stack with no other files in the folder
    %       3. Load tiffstack into matlab as a PIVdroplets object
    %       4. Run the find_pos_hough code to find droplet positions
    %       5. Run the track_obj code to connect droplet positions
    %       6. Iterate step 4 & 5 to optimize parameters
    %       6b. Check trace quality by running the code plot_tracks
    %       7. Run the code smooth_trace to connect any missing links
    %       8. Run the code output_PIV_mask to do just that. A file named
    %       contour will be saved in the current folder
    %       9. Open the PIVlab application and select the folder with
    %       individual tiff files to load in
    %       10. Load data separated as 1-2, 3-4. Load the contour file by
    %       using "load contour", NOT "load external contour"
    %       11. Set the analysis parameters as you prefer them
    %       12. Run the PIVlab file, postprocess & export as .mat file
    %       13. Close PIVlab & open the PIVlab mat file you exported in the
    %       workspace. Delete all files except x,y,u_original and
    %       v_original
    %       14. Load those files into the PIVdroplet object by running the
    %       code add_PIVdata_new
    %       15. Run code average_profile_rot2(obj). This takes a while and
    %       will average all flow profiles while translating and rotating
    %       the center of the droplet such that they overlap in each frame
    %       and the swimming direction is the same.
    %       16. Plot the profile with streamlines using the
    %       plot_streamlines_average_rot function
    %
    properties
        Flows %1xNOF cell structure that contains a profile_class object in each cell
        average_flow %1x1 structure of type profileclass containing the average flow profile
        average_flow_trans %1x1 structure of type profileclass containing the average flow profile after only translation
        average_flow_rot %1x1 structure of type profileclass containing the average flow profile after translation and rotation
        average_flow_rotmintrans %1x1 structure of type profileclass containing the average rotated profile with the average translated profile substracted
        trace_short %matches the size of the flow profiles and ignores every other frame
        vellistrot %list of xvel and yvel of droplets rotated in each frame to align with x axis
    end
    methods
        function obj=PIVdroplets(fname,parameters,identity)
            %Initializes object identical to droplets
            
            if nargin==1
                parameters=[];
            elseif nargin<1
                parameters=[];
                fname=[];
            end
            if isempty(parameters)
                initialize=1;
            else
                initialize=0;
            end
            if nargin<3
                identity='y';
            end
            
            obj=obj@droplets(fname,parameters,identity);
        end
        function obj=add_PIVdata(obj)
            %Takes in n files of .mat type containing PIVlab output, saves
            %them as profile_class and stores them into the Flows property
            %of this object.
            %assumes filename of .mat file has type "name_####.mat" where
            %#### is a 4 digit number starting at 0001.
            
            [filen, foldername_temp] = uigetfile('*.mat','Select first .mat file with ID 0001');
            filen_sub=filen(1:length(filen)-8);
            
            if isempty(obj.pos)
                obj=find_pos_hough(obj);
            end
            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            %Create an empty cell structure for storing flow
            flow_temp=cell(obj.NOF,1);
            
            for k=1:obj.NOF
                %Find name of file with number k
                filename2=[foldername_temp filen_sub num2str(k,'%04d') '.mat'];
                check=1;
                try filetemp=open(filename2);
                catch
                    check=0;
                    disp(['A file with filename "' filename2 '" could not be opened.'])
                end
                
                %Open the x y u and v vector in object contain in mat file
                %number k
                if check==1
                    u=filetemp.u;
                    v=filetemp.v;
                    x=filetemp.x;
                    y=filetemp.y;
                    
                    %Create an object of type profile_class and store in flow
                    flow_temp{k,1}=profile_class6(x,y,u,v);
                else
                    flow_temp{k,1}=NaN;
                end
            end %End forloop
            obj.Flows=flow_temp;
        end
        function obj=contract_trace(obj)
            obj.trace_short=[];
            for p=1:obj.NOP
                tr_p=obj.tr(find(obj.tr(:,4)==p),1:4);
                tr_p=tr_p(1:2:end-1,:);
                obj.trace_short=[obj.trace_short ; tr_p];
            end
        end
        function obj=add_PIVdata_new(obj,xlist,ylist,ulist,vlist)
            flow_temp=cell(floor(obj.NOF/2),1);
            for i=1:floor(obj.NOF/2)
                x=xlist{i};
                y=ylist{i};
                u=ulist{i};
                v=vlist{i};
                flow_temp{i,1}=profile_class6(x,y,u,v);
            end
            obj.Flows=flow_temp;
        end
        function output=export_average_profile_rot(obj)
            %Exports the average profile (rot) in a format that pivlab can
            %import so that streamlines can be drawn
            
            begincount=floor((length(obj.average_flow_rot.x)-63)/2);
            output{1,1}=obj.average_flow_rot.x(begincount:begincount+62,begincount:begincount+62);
            output{1,1}=output{1,1}-min(min(output{1,1}))+1;
            output{2,1}=obj.average_flow_rot.y(begincount:begincount+62,begincount:begincount+62);
            output{2,1}=output{2,1}-min(min(output{2,1}))+1;
            output{3,1}=obj.average_flow_rot.u(begincount:begincount+62,begincount:begincount+62);
            output{4,1}=obj.average_flow_rot.v(begincount:begincount+62,begincount:begincount+62);
            output{5,1}=ones(size(output{1,1}));
            output{6,1}=[];
        end
        function obj=average_profile(obj,selection)
            %Makes a list of profiles
            %Averages them, and returns the average flow profile
            %selection must be a list of numbers indicating which frames to
            %consider. eg [1:100].
            
            if nargin<2
                selection=0;
            end
            
            %check if flow information is present
            if isempty(obj.Flows)
                obj=add_PIVdata(obj);
            end
            
            %make a list of all profiles in the stack;
            listofprofiles=obj.Flows{1,1};
            for n=2:floor(obj.NOF/2)
                new_profile=obj.Flows{n,1};
                try listofprofiles=[listofprofiles new_profile];
                catch
                end
            end
            
            if selection~=0
                %take a selection from the list
                if max(selection)>length(listofprofiles)
                    error('Selection is out of range, please provide a valid selection.')
                end
                try listofprofiles=listofprofiles(selection);
                catch
                    error('Your selection does not have the right format or indicates values that are out of range. Please provide a valid selection.')
                end
            end
            
            %average the profiles
            averageflow=average_profiles_fast(listofprofiles);
            obj.average_flow=averageflow;
            
            %plot the obtained average;
            plot_flow(obj);
            
        end
        function obj=average_profile_rot(obj,selection)
            %Takes the flow profiles around a swimming droplet, rotates and
            %translates such that all droplets sit on the origin and move
            %in the same direction
            %selection is optional and can be a list of frames of interest
            %organized as [start:end] or [start middle end];
            
            if nargin<2
                selection=0;
            end
            
            %Make sure tracking data are presesnt
            if isempty(obj.trace_short)
                obj=contract_trace(obj);
            end
            
            %Make sure flows are present;
            if isempty(obj.Flows)
                obj=add_PIVdata(obj);
            end
            
            % Get particle trace
            %trace=obj.trace_short;
            trace=obj.tr;
            
            %Check if selection option is used
            if selection==0
                nframes=floor(obj.NOF/2);
            else
                nframes=max(selection);
            end
            
            %define u2,v2,X,Y (velocity and coordinate components particle)
            U=NaN(nframes,1);
            V=NaN(nframes,1);
            X=NaN(nframes,1);
            Y=NaN(nframes,1);
            Uactual=NaN(nframes,1);
            Vactual=NaN(nframes,1);
            
            % Get particle location and velocity
            for k=2:nframes
                if any(selection==k) | selection==0
                    %calculate x component motion droplet
                    k+nframes;
                    U(k,:)=trace(k+nframes,1)-trace(k,1);%actually the distance between particle 1 and 2
                    Uactual(k,:)=(trace(k,1)-trace(k-1,1)+trace(nframes+k,1)-trace(nframes+k-1,1))/2;
                    Vactual(k,:)=(trace(k,2)-trace(k-1,2)+trace(nframes+k,2)-trace(nframes+k-1,2))/2;
                    %calculate y component motion droplet
                    V(k,:)=trace(k+nframes,2)-trace(k,2);
                    X(k,:)=(trace(k,1)+trace(k+nframes,1))/2;
                    Y(k,:)=(trace(k,2)+trace(k+nframes,2))/2;
                end
            end
            
            % Create list of translated and rotated profiles
            profile_list=[];
            velocitylist=NaN(nframes,2);
            for k=2:nframes
                %collect matfile of framenumber k
                framek=obj.Flows{k,1};
                
                %Calculate rotation angle
                %Note: tan alpha is y-component of motion divided by x component of motion
                %Minus sign because we want to rotate back to the x-axis.
                if U(k)>0
                    ang = -atan(V(k)/U(k));
                else
                    ang = pi - atan(V(k)/U(k));
                end
                
                framek=translate_profile(framek,-X(k),-Y(k));
                framek=rotate_profile(framek,ang);
                rotvec=[cos(ang) -sin(ang) ;  sin(ang) cos(ang)]*[Uactual(k) ; Vactual(k)];
                Urot=rotvec(1);
                Vrot=rotvec(2);
                velocitylist(k,1:2)=[Urot Vrot];
                profile_list=[profile_list framek];
            end
            
            obj.vellistrot=velocitylist;
            meanU=nanmean(velocitylist(:,1));
            meanV=nanmean(velocitylist(:,2));
            
            % Average and plot average
            I=obj.im;
            averaged_profile=average_profiles_fast(profile_list);
            figure
            hold on
            box on
            colormap gray
            imagesc(I);
            plot_profile(averaged_profile)
            hold off
            
            %Save average into object
            obj.average_flow_rot=averaged_profile;
            
            
        end
        function obj=average_profile_rot2(obj,selection)
            %Takes the flow profiles around a swimming droplet, rotates and
            %translates such that all droplets sit on the origin and such
            %that THE VECTOR GOING FROM DROPLET 1 TO DROPLET 2 IS THE SAME AS IN FRAME 1
            %This is different from average_profile_rot (not 2) in that the
            %speed of the cluster is not used, but rather it's orientation
            %selection is optional and can be a list of frames of interest
            %organized as [start:end] or [start middle end];
            if nargin<2
                selection=0;
            end
            %Make sure tracking data are presesnt
            if isempty(obj.trace_short)
                obj=contract_trace(obj);
            end
            %Make sure flows are present;
            if isempty(obj.Flows)
                obj=add_PIVdata_new(obj);
            end
            % Get particle trace
            %trace=obj.trace_short;
            trace=obj.tr;
            %Check if selection option is used
            if selection==0
                nframes=floor(obj.NOF/2);
            else
                nframes=max(selection);
            end
            %define u2,v2,X,Y (direction and center of the particles)
            U=NaN(nframes,1);
            V=NaN(nframes,1);
            Xmean=NaN(nframes,1);
            Ymean=NaN(nframes,1);
            for k=1:nframes
            %for k=1:2
                if any(selection==k) | selection==0
                    
                    %calculate the x and y component of the cluster
                    %position
                    X=NaN(1,obj.NOP);
                    Y=NaN(1,obj.NOP);
                    for p=1:obj.NOP
                        %get the x and y position of each particle in the
                        %cluster
                        tracep=obj.tr(find(obj.tr(:,4)==p),1:2);
                        %XYp(p,1:2)=trace(k+(p-1)*nframes,1:2);
                        X(1,p)=tracep(k,1);
                        Y(1,p)=tracep(k,2);
                    end
                    Xmean(k,1)=nanmean(X(1,:));
                    Ymean(k,1)=nanmean(Y(1,:));
                    
                    %If cluster consists of 2 particles, you can use their
                    %relative position to reorient the cluster instead of
                    %the velocity
                    if obj.NOP==2
                        U(:,1)=X(1,2)-X(1,1);
                        V(:,1)=Y(1,2)-Y(1,1);
                    %else you have to use the motion of the droplet to
                    %define its direction
                    else
                        U(k,:)=trace(k+nframes,1)-trace(k,1);%actually the x distance between particle 1 and 2
                        %calculate y component motion droplet
                        V(k,:)=trace(k+nframes,2)-trace(k,2);
                    end
                end
            end
            disp(['Analyzing ' num2str(nframes) ' nframes.']);
            % Create list of translated and rotated profiles
            profile_list=[];
            for k=1:nframes
            %for k=1:2
                %collect matfile of framenumber k
                framek=obj.Flows{k,1};
                
                %Calculate rotation angle
                %Minus sign because we want to rotate back to the x-axis.
                ang=-atan2(V(k),U(k));
                
                %Translate the profile so that the center of mass is always
                %on [0,0]
                framek=translate_profile(framek,-Xmean(k,1),-Ymean(k,1));
                %Rotate the profile so that the vector of particle 1 -
                %particle 2 is along the positive x axis
                framek=rotate_profile(framek,ang);
                %plot_profile(framek)
                %save the rotated translated profile in a list
                profile_list=[profile_list framek];
            end
            
            %Average and plot average
            I=obj.im;
            averaged_profile=average_profiles_fast(profile_list);
            figure
            hold on
            box on
            colormap gray
            imagesc(I);
            plot_profile(averaged_profile)
            hold off
            
            %Save average into object
            obj.average_flow_rot=averaged_profile;
            
            
        end
        function obj=average_profile_trans(obj,selection)
            %Takes the flow profiles around a swimming droplet, rotates and
            %translates such that all droplets sit on the origin and move
            %in the same direction
            %selection is optional and can be a list of frames of interest
            %organized as [start:end] or [start middle end];
            
            if nargin<2
                selection=0;
            end
            
            %Make sure tracking data are presesnt
            if isempty(obj.trace_short)
                obj=contract_trace(obj);
            end
            
            %Make sure flows are present;
            if isempty(obj.Flows)
                obj=add_PIVdata(obj);
            end
            
            % Get particle trace
            %trace=obj.trace_short;
            trace=obj.tr;
            
            %Check if selection option is used
            if selection==0
                nframes=floor(obj.NOF/2);
            else
                nframes=max(selection);
            end
            
            %define u2,v2,X,Y (velocity and coordinate components particle)
            U=NaN(nframes,1);
            V=NaN(nframes,1);
            X=NaN(nframes,1);
            Y=NaN(nframes,1);
            
            % Get particle location and velocity
            for k=1:nframes
                if any(selection==k) | selection==0
                    %calculate x component motion droplet
                    U(k,:)=trace(k+1,1)-trace(k,1);
                    %calculate y component motion droplet
                    V(k,:)=trace(k+1,2)-trace(k,2);
                    X(k,:)=(trace(k,1)+trace(k+nframes,1))/2;
                    Y(k,:)=(trace(k,2)+trace(k+nframes,2))/2;
                end
            end
            
            % Create list of translated and rotated profiles
            profile_list=[];
            for k=1:nframes
                %collect matfile of framenumber k
                framek=obj.Flows{k,1};
                framek=translate_profile(framek,-X(k),-Y(k));
                profile_list=[profile_list framek];
            end
            
            % Average and plot average
            averaged_profile=average_profiles_fast(profile_list);
            plot_profile(averaged_profile)
            
            %Save average into object
            obj.average_flow_trans=averaged_profile;
        end
        function obj=substract_rot_trans(obj)
            %Takes the rotated average profile and substracts the
            %translated average profile.
            
            %make sure both profiles exist, otherwise fill them
            if isempty(obj.average_flow_rot)
                obj=average_profile_rot(obj);
            end
            if isempty(obj.average_flow_trans)
                obj=average_profile_trans(obj);
            end
            
            %get the two profiles to be substracted
            rot_temp=obj.average_flow_rot;
            trans_temp=obj.average_flow_trans;
            
            %get a new grid that fits both
            new_grid=get_new_grid([rot_temp,trans_temp]);
            
            %fit both profiles on that grid
            rot_temp2=reform_grid(rot_temp,new_grid.x,new_grid.y);
            trans_temp2=reform_grid(trans_temp,new_grid.x,new_grid.y);
            
            %substract the profiles and plot
            rot_substract=substract_profile(rot_temp2,trans_temp2);
            plot_profile(rot_substract,[-1000 1000 -1000 1000]);
            
            %save profile in object
            obj.average_flow_rotmintrans=rot_substract;
            
        end
        function plot_streamlines_average_rot(obj,npoints,r)
            %plots a figure with streamlines around the average
            %Collect data
            limval=300;
            close all
            
            x=obj.average_flow_rot.x;
            y=obj.average_flow_rot.y;
            u=obj.average_flow_rot.u*100;
            v=obj.average_flow_rot.v*100;
            
            %define startpoints for streamlines around the central object
            %in radius r
            startx=NaN(1,npoints);
            starty=NaN(1,npoints);
            for j=1:npoints
                startx(1,j)=cos(2*pi*j/npoints)*r;
                starty(1,j)=sin(2*pi*j/npoints)*r;
            end
            
            %size(startx)
            %startx=[startx' ;150 ;150 ; 150 ; 150;150; 150];
            %starty=[starty' ; -250; -200 ;-150;150 ;200; 250];
            %startx=startx';
            %starty=starty';
            startx=[startx -250 -250 -250 ];
            starty=[starty 0 50 -50];
            
             %define a second layer of startpoints for streamlines further around the central object
            %in radius 3r
            startx2=NaN(1,npoints);
            starty2=NaN(1,npoints);
            for j=1:npoints*3
                startx2(1,j)=cos(2*pi*(j+15)/npoints)*3*r;
                starty2(1,j)=sin(2*pi*(j+15)/npoints)*3*r;
            end
            
            figure(1)
            hold on
            box on
            streamline(x,y,u,v,startx,starty)
            streamline(x,y,-u,-v,startx,starty)
            %streamline(x,y,u,v,startx2,starty2)
            %streamline(x,y,-u,-v,startx2,starty2)
            pbaspect([1 1 1])
            xlim([-limval limval]);
            ylim([-limval limval]);
            hold off
            fig1=figure(1);
            fig1.Renderer='Painters';
            
            figure(2)
            hold on
            box on
            quiver(x,y,u,v,2,'LineWidth',1,'Color','r')
            pbaspect([1 1 1])
            xlim([-limval limval]);
            ylim([-limval limval]);
            hold off
            fig2=figure(2);
            fig2.Renderer='Painters';
        end
        function plot_streamlines_average(obj,npoints,r)
            %plots a figure with streamlines around the average
            %Collect data
            limval=300;
            close all
            
            x=obj.average_flow.x;
            y=obj.average_flow.y;
            u=obj.average_flow.u*100;
            v=obj.average_flow.v*100;
            
            %define startpoints for streamlines around the central object
            %in radius r
            startx=NaN(1,npoints);
            starty=NaN(1,npoints);
            for j=1:npoints
                startx(1,j)=cos(2*pi*j/npoints)*r;
                starty(1,j)=sin(2*pi*j/npoints)*r;
            end
            
            %size(startx)
            %startx=[startx' ;150 ;150 ; 150 ; 150;150; 150];
            %starty=[starty' ; -250; -200 ;-150;150 ;200; 250];
            %startx=startx';
            %starty=starty';
            startx=[startx -250 -250 -250 ];
            starty=[starty 0 50 -50];
            
             %define a second layer of startpoints for streamlines further around the central object
            %in radius 3r
            startx2=NaN(1,npoints);
            starty2=NaN(1,npoints);
            for j=1:npoints*3
                startx2(1,j)=cos(2*pi*(j+15)/npoints)*3*r;
                starty2(1,j)=sin(2*pi*(j+15)/npoints)*3*r;
            end
            
            figure(1)
            hold on
            box on
            streamline(x,y,u,v,startx,starty)
            streamline(x,y,-u,-v,startx,starty)
            %streamline(x,y,u,v,startx2,starty2)
            %streamline(x,y,-u,-v,startx2,starty2)
            pbaspect([1 1 1])
            xlim([-limval limval]);
            ylim([-limval limval]);
            hold off
            fig1=figure(1);
            fig1.Renderer='Painters';
            
            figure(2)
            hold on
            box on
            quiver(x,y,u,v,2,'LineWidth',1,'Color','r')
            pbaspect([1 1 1])
            xlim([-limval limval]);
            ylim([-limval limval]);
            hold off
            fig2=figure(2);
            fig2.Renderer='Painters';
        end
        function output=plot_streamlines_magnitudecolored_rot(obj,npoints,r)
            %plots a figure with streamlines around the average
            %Collect data
            limval=300;
            close all
            
            x=obj.average_flow_rot.x;
            y=obj.average_flow_rot.y;
            u=obj.average_flow_rot.u;
            v=obj.average_flow_rot.v;
            
            I=sqrt(v.^2+u.^2);
            
            %define startpoints for streamlines around the central object
            %in radius r
            startx=NaN(1,npoints);
            starty=NaN(1,npoints);
            for j=1:npoints
                startx(1,j)=cos(2*pi*j/npoints)*r;
                starty(1,j)=sin(2*pi*j/npoints)*r;
            end
            
            %size(startx)
            %startx=[startx' ;150 ;150 ; 150 ; 150;150; 150];
            %starty=[starty' ; -250; -200 ;-150;150 ;200; 250];
            %startx=startx';
            %starty=starty';
            %startx=[startx -250 -250 -250 ];
            %starty=[starty 0 50 -50];
            
             %define a second layer of startpoints for streamlines further around the central object
            %in radius 3r
            startx2=NaN(1,npoints);
            starty2=NaN(1,npoints);
            for j=1:npoints*3
                startx2(1,j)=cos(2*pi*(j+15)/npoints)*3*r;
                starty2(1,j)=sin(2*pi*(j+15)/npoints)*3*r;
            end
            
            figure(1)
            hold on
            box on
            streamline(x,y,u,v,startx,starty)
            streamline(x,y,-u,-v,startx,starty)
            %streamline(x,y,u,v,startx2,starty2)
            %streamline(x,y,-u,-v,startx2,starty2)
            pbaspect([1 1 1])
            xlim([-limval limval]);
            ylim([-limval limval]);
            hold off
            fig1=figure(1);
            fig1.Renderer='Painters';
            
            figure(2)
            hold on
            box on
            quiverC2D(x,y,u,v,2,'LineWidth',1.5,'MaxHeadSize',1)
            pbaspect([1 1 1])
            xlim([-limval limval]);
            ylim([-limval limval]);
            hold off
            fig2=figure(2);
            fig2.Renderer='Painters';
            
            trace=obj.trace_short;
            nframes=floor(obj.NOF/2);
            U=trace(1+nframes,1)-trace(1,1);
            V=trace(1+nframes,2)-trace(1,2);
            ang=-atan2(V,U);
            im1=imrotate(obj.im,ang);
            %figure(3)
            %hold on
            %box on
            %imagesc(-trace(1,1),-trace(2,2),im1);
            %streamline(x,y,u,v,startx,starty)
            %streamline(x,y,-u,-v,startx,starty)
            %quiverC2D(x,y,u,v,2,'LineWidth',1)
            %pbaspect([1 1 1])
            %xlim([-limval limval]);
            %ylim([-limval limval]);
            %hold off
            %fig3=figure(3);
            %fig3.Renderer='Painters';
            output=I;
        end
        function output=plot_streamlines_magnitudecolored(obj,npoints,r)
            %plots a figure with streamlines around the average
            %Collect data
            limval=250;
            close all
            
            x=obj.average_flow.x;
            y=obj.average_flow.y;
            u=obj.average_flow.u;
            v=obj.average_flow.v;
            
            I=sqrt(v.^2+u.^2);
            
            %define startpoints for streamlines around the central object
            %in radius r
            startx=NaN(1,npoints);
            starty=NaN(1,npoints);
            for j=1:npoints
                startx(1,j)=cos(2*pi*j/npoints)*r;
                starty(1,j)=sin(2*pi*j/npoints)*r;
            end
            
            %size(startx)
            %startx=[startx' ;150 ;150 ; 150 ; 150;150; 150];
            %starty=[starty' ; -250; -200 ;-150;150 ;200; 250];
            %startx=startx';
            %starty=starty';
            %startx=[startx -250 -250 -250 ];
            %starty=[starty 0 50 -50];
            
             %define a second layer of startpoints for streamlines further around the central object
            %in radius 3r
            startx2=NaN(1,npoints);
            starty2=NaN(1,npoints);
            for j=1:npoints*3
                startx2(1,j)=cos(2*pi*(j+15)/npoints)*3*r;
                starty2(1,j)=sin(2*pi*(j+15)/npoints)*3*r;
            end
            
            figure(1)
            hold on
            box on
            streamline(x,y,u,v,startx,starty)
            streamline(x,y,-u,-v,startx,starty)
            %streamline(x,y,u,v,startx2,starty2)
            %streamline(x,y,-u,-v,startx2,starty2)
            pbaspect([1 1 1])
            xlim([-limval limval]);
            ylim([-limval limval]);
            hold off
            fig1=figure(1);
            fig1.Renderer='Painters';
            
            figure(2)
            hold on
            box on
            quiverC2D(x,y,u,v,2,'LineWidth',1)
            pbaspect([1 1 1])
            xlim([-limval limval]);
            ylim([-limval limval]);
            hold off
            fig2=figure(2);
            fig2.Renderer='Painters';
            
            trace=obj.tr;
            nframes=floor(obj.NOF/2);
            U=trace(1+nframes,1)-trace(1,1);
            V=trace(1+nframes,2)-trace(1,2);
            ang=-atan2(V,U);
            im1=imrotate(obj.im,ang);
            %figure(3)
            %hold on
            %box on
            %imagesc(-trace(1,1),-trace(2,2),im1);
            %streamline(x,y,u,v,startx,starty)
            %streamline(x,y,-u,-v,startx,starty)
            %quiverC2D(x,y,u,v,2,'LineWidth',1)
            %pbaspect([1 1 1])
            %xlim([-limval limval]);
            %ylim([-limval limval]);
            %hold off
            %fig3=figure(3);
            %fig3.Renderer='Painters';
            output=I;
        end
        function obj=plot_flow(obj,imcheck)
            %plots the average flow profile
            %imcheck controls whether (1) or not (0) image is shown under flow
            %default is not (0);
            
            if nargin<2
                imcheck=0;
            end
            
            %Check if flow profile information is present
            if isempty(obj.average_flow)
                obj=average_profile(obj);
            end
            
            %Unpack flow profile information
            averageflow=obj.average_flow;
            %rotate to have top be top
            averageflow=rotate_profile(averageflow,pi);
            %translate so that image starts at 0 0 in bottom corner
            imsize=size(obj.im);
            averageflow=translate_profile(averageflow,imsize(2),imsize(1));
            %unpack object
            x=averageflow.x;
            y=averageflow.y;
            u=averageflow.u;
            v=averageflow.v;
            
            %Flip image over y
            image=rot90(obj.im,2);
            %image=flipdim(obj.im,2);
            
            %Assert image is 8-bit
            if length(size(image))>2
                disp('Image is not of type 8-bit. Please convert to 8-bit first')
                %Convert to grayscale
                image = rgb2gray(image);
            end
            
            %Increase image brightness
            image=image+(255-max(max(image)));
            
            
            %plot object
            figure
            hold on
            if imcheck~=0
                colormap('gray'), imagesc(image);
            end
            box on
            axis([0 max(max(x)) 0 max(max(y))]);
            pbaspect([max(max(x))/max(max(y)) 1 1]);
            quiver(x,y,u,v,1,'b');
            hold off
            
        end
        function obj=plot_flow_rot(obj,imcheck)
            %plots the average flow profile
            %imcheck controls whether (1) or not (0) image is shown under flow
            %default is not (0);
            
            if nargin<2
                imcheck=0;
            end
            
            %Check if flow profile information is present
            if isempty(obj.average_flow_rot)
                obj=average_profile_rot(obj);
            end
            
            %Unpack flow profile information
            averageflow=obj.average_flow_rot;
            %rotate to have top be top
            averageflow=rotate_profile(averageflow,pi);
            %translate so that image starts at 0 0 in bottom corner
            imsize=size(obj.im);
            %averageflow=translate_profile(averageflow,imsize(2),imsize(1));
            %unpack object
            x=averageflow.x;
            y=averageflow.y;
            u=averageflow.u;
            v=averageflow.v;
            
            %Flip image over y
            %image=rot90(obj.im,2);
            image=flipdim(obj.im,1);
            
            %Assert image is 8-bit
            if length(size(image))>2
                disp('Image is not of type 8-bit. Please convert to 8-bit first')
                %Convert to grayscale
                image = rgb2gray(image);
            end
            
            %Increase image brightness
            image=image+(255-max(max(image)));
            
            velocitylist=obj.vellistrot;
            up1=nanmean(velocitylist(:,1))
            vp1=nanmean(velocitylist(:,2))
            sizeu=size(u);
            sizev=size(v);
            up=zeros(sizeu);
            vp=zeros(sizev);
            up(round(sizeu(1)/2),round(sizeu(2)/2))=up1;
            vp(round(sizev(1)/2),round(sizev(2)/2))=vp1;
            
            scaleval=20;
            %plot object
            figure
            hold on
            box on
            if imcheck~=0
                colormap('gray'), imshow(image);
            end
            %axis([0 max(max(x)) 0 max(max(y))]);
            %pbaspect([max(max(x))/max(max(y)) 1 1]);
            quiver(x,y,scaleval*u,scaleval*v,0,'b');
            quiver(x,y,scaleval*up,scaleval*vp,0,'r');
            hold off
            
        end
        function obj=plot_flow_trans(obj,imcheck)
            %plots the average flow profile
            %imcheck controls whether (1) or not (0) image is shown under flow
            %default is not (0);
            
            if nargin<2
                imcheck=0;
            end
            
            %Check if flow profile information is present
            if isempty(obj.average_flow_trans)
                obj=average_profile_trans(obj);
            end
            
            %Unpack flow profile information
            averageflow=obj.average_flow_trans;
            %rotate to have top be top
            averageflow=rotate_profile(averageflow,pi);
            %translate so that image starts at 0 0 in bottom corner
            imsize=size(obj.im);
            averageflow=translate_profile(averageflow,imsize(2),imsize(1));
            %unpack object
            x=averageflow.x;
            y=averageflow.y;
            u=averageflow.u;
            v=averageflow.v;
            
            %Flip image over y
            image=rot90(obj.im,2);
            
            %Assert image is 8-bit
            if length(size(image))>2
                disp('Image is not of type 8-bit. Please convert to 8-bit first')
                %Convert to grayscale
                image = rgb2gray(image);
            end
            
            %Increase image brightness
            image=image+(255-max(max(image)));
            
            %plot object
            figure
            hold on
            if imcheck~=0
                colormap('gray'), imagesc(image);
            end
            box on
            %axis([0 max(max(x)) 0 max(max(y))]);
            axis([310 1280 80 580])
            pbaspect([(1280-310)/(500) 1 1])
            %pbaspect([max(max(x))/max(max(y)) 1 1]);
            quiver(x,y,u,v,2,'b');
            hold off
            
        end
        function output=plot_velocity_contour(obj)
            %plots a contourmap of velocity around droplets
            %requires either average_flow or average_flow_trans as input

            %close figures
            close all
            %indicate value for cropping data
            startval=6;
            
            
            %Unpack flow profile information
            if isempty(obj.average_flow)    
                averageflow=obj.average_flow_trans;
            else
                averageflow=obj.average_flow;
            end
            %rotate to have top be top
            averageflow=rotate_profile(averageflow,pi);
            %crop data
            x=averageflow.x(1+startval:end-startval,1+startval:end-startval);
            y=averageflow.y(1+startval:end-startval,1+startval:end-startval);
            u=averageflow.u(1+startval:end-startval,1+startval:end-startval);
            v=averageflow.v(1+startval:end-startval,1+startval:end-startval);
            
            %calculate speed
            I=sqrt(v.^2+u.^2)*obj.scale*obj.framerate;
            
            %plot image
            figure(1)
            hold on
            box on            
            contourf(x,y,I)
            colorbar('southoutside');
            %quiverC2D(x,y,u,v,2,'LineWidth',1)
            pbaspect([1 1 1])
            xlim([-400 400]);
            ylim([-200 600]);
            hold off
            fig1=figure(1);
            fig1.Renderer='Painters';
            output=I;
        end
        function obj=substract_background_flow(obj,value)
            %Value is the magnitude of the background flow of form [x y];
            
            if isempty(obj.average_flow)
                obj=average_profile(obj);
            end
            averageflow=obj.average_flow;
            
            new_grid.x=averageflow.x;
            new_grid.y=averageflow.y;
            new_grid.u=averageflow.x-averageflow.x+(1*value(1));
            new_grid.v=averageflow.x-averageflow.x+(1*value(2));
            
            averageflow=substract_profile(averageflow,new_grid);
            obj.average_flow=averageflow;
            plot_flow(obj);
        end
        function plot_polar_decay(obj,frame)
            %plots the radial component of the speed as a function of
            %distance to droplet center for framenumber 'frame'
            %if framenumber is not provided, uses the average
            %if framenumber is not provided and average is not calculated,
            %uses the first frame
            %only works if only one droplet is present
            
            if obj.NOP>1
                error('This function can only handle one droplet at a time')
            end
            
            %make sure flows are provided
            if isempty(obj.Flows)
                obj=add_PIVdata(obj);
            end
            
            %make sure trace is provided
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            %get the profile of interest
            if nargin==1
                if isempty(obj.average_flow_rot)
                    framen=obj.Flows{1};
                    centerx=obj.tr(1,1);
                    centery=obj.tr(1,2);
                else
                    framen=obj.average_flow_rot;
                    centerx=0;
                    centery=0;
                end
            else
                framen=obj.Flows{1};
                centerx=obj.tr(1,1);
                centery=obj.tr(1,2);
            end
            
            
            %get decay information [distance v_r v_theta] using
            %profile_class6
            decay=get_radial_decay_polar(framen,centerx,centery);
            
            figure
            hold on
            box on
            title('radial component');
            xlabel('r (\mum)','FontSize',18)
            ylabel('v_r (\mum/s)','FontSize',18)
            plot(decay(:,1)*obj.scale,decay(:,2)*obj.scale*obj.framerate,'ko');
            hold off
            
            figure
            hold on
            box on
            title('angular component');
            xlabel('r (\mum)','FontSize',18)
            ylabel('v_\theta (\mum/s)','FontSize',18)
            plot(decay(:,1)*obj.scale,decay(:,3)*obj.scale*obj.framerate,'ko');
            hold off
            
            
        end
        function plot_trace_on_movie_flow(obj,start_frame, end_frame)
            %A is the total stack of images corresponding to the track
            %Displays a movie of the droplets in each frame with the
            %corresponding position on top.
            %Settings for video creation
            %Saves the video to the indicated foldername
            %Uses a framerate of 7;
            
            if isnan(obj.im)
                error('There is no imagestack loaded in the object, so no trajectory can be plotted')
            end
            %if trace does not exist, trace particles first
            if isempty(obj.trace_short)
                obj=contract_trace(obj);
            end
            
            %if flows do not exist, add PIV data first
            if isempty(obj.Flows)
                obj=add_PIVdata(obj);
            end
            
            linesize=2;
            %speed prefactor
            factor=1;
            %speed linewidth
            linewidth=1;
            nframes=floor(obj.NOF/2);
            
            if nargin<3
                end_frame=nframes-1;
            end
            if nargin<2
                end_frame=nframes-1;
                start_frame=1;
            end
            
            if end_frame>nframes
                end_frame=nframes-1;
                disp('The last frame you indicated is larger than the number of frames. We are using the whole movie.')
            end
            
            nparts=obj.NOP;
            trace=obj.trace_short;
            fname=[obj.foldername obj.filename '.tif'];
            
            %define v2: list of y compoments of particle motion
            %define u2: list of x components of particle motion
            u2=NaN(nframes,nparts);
            v2=NaN(nframes,nparts);
            
            %Get a quiverplot of the direction of motion of each particle
            for i=1:nparts
                trj=trace(find(trace(:,4)==i),1:3);
                for k=1:nframes-1
                    x2(k,i)=trj(k,1);
                    y2(k,i)=trj(k,2);
                    %calculate x component motion particle
                    u2(k,i)=trj(k+1,1)-trj(k,1);
                    %calculate y component motion particle
                    v2(k,i)=trj(k+1,2)-trj(k,2);
                end
            end
            
            qqmax=ceil(end_frame/50);
            moviefolder=obj.foldername;
            
            for qq=start_frame:qqmax
                
                moviename=[obj.filename '_trackmovie_' num2str(qq) '.avi'];
                
                t=(qq-1)*50+1;
                
                I=imread(fname,t);
                %define size image for aspect ratio
                sizeI=size(I);
                
                
                %get flow data of frame 1
                X=obj.Flows{t,1}.x;
                Y=obj.Flows{t,1}.y;
                U=obj.Flows{t,1}.u;
                V=obj.Flows{t,1}.v;
                
                % open video
                name_video=[ moviefolder '\' moviename];
                video = VideoWriter(name_video,'Uncompressed AVI');
                video.FrameRate = 7;
                open(video);
                
                %Plot new figure of frame 1
                figure
                hold on
                %show original image in grayscale
                colormap gray
                imagesc(I);
                %show vectorplot of particle motion in red
                for i=1:nparts
                    quiver(x2(t,i),y2(t,i),factor*u2(t,i),factor*v2(t,i),linesize,'Color','r','LineWidth',linewidth)
                end
                %show vectorplot of liquid flow in blue
                quiver(X,Y,U,V,1,'Color','b','LineWidth',linewidth)
                
                %Settings aspect ratios
                set(gca,'Ylim',[0 sizeI(1)],'Xlim', [0 sizeI(2)])
                set(gca,'Ytick',[0:200:sizeI(2)],'Xtick',[0:300:sizeI(1)])
                pbaspect([sizeI(1)/sizeI(1) sizeI(1)/sizeI(2) 1]);
                set(gca,'YDir','normal');
                
                %Settings axis lables
                xlabel('pixels','Fontname','Helvetica','Fontsize',14);
                ylabel('pixels','Fontname','Helvetica','Fontsize',14);
                
                %Settings frame formatting
                set(gca,'Box','on','XMinorTick','on','Yminortick','on');
                set(gca,'TickLength',[0.025 0.035]);
                set(gca,'LineWidth',1.5);
                set(gcf,'Color',[1 1 1]);
                set(gcf, 'InvertHardCopy', 'off');
                set(gca,'layer','top');
                set(gcf,'Paperunits','inches','Papersize',[6.5 6],'Paperposition',[0 0 6.5 6]);
                
                
                %next plot use same layout as this plot (needed for video)
                set(gca,'nextplot','replacechildren');
                %for all frames in folder
                tmax=t+50;
                if tmax>nframes && qq~=1
                    tmax=nframes;
                end
                if tmax>end_frame && qq==1
                    tmax=end_frame;
                end
                
                for k=t+1:tmax-1
                    %read image
                    fname=[obj.foldername obj.filename '.tif'];
                    Ik=imread(fname,2*k);
                    
                    %read flow profiles
                    X=obj.Flows{k,1}.x;
                    Y=obj.Flows{k,1}.y;
                    U=obj.Flows{k,1}.u;
                    V=obj.Flows{k,1}.v;
                    
                    % Print number of frames still to process
                    if ~mod((nframes-k)/100,1)
                        disp(['Still to process: ' num2str(nframes - k,'%03d') ' frames'])
                    end
                    %within image
                    hold on
                    %plot in grayscale
                    colormap gray
                    %show tif image
                    imagesc(Ik)
                    %plot direction of particle
                    for i=1:nparts
                        quiver(x2(k,i),y2(k,i),factor*u2(k,i),factor*v2(k,i),linesize,'Color','r','LineWidth',linewidth)
                    end
                    %plot flow profile
                    quiver(X,Y,U,V,1,'Color','b','LineWidth',linewidth)
                    
                    %make sure formatting in move is same as above
                    set(gca,'nextplot','replacechildren');
                    
                    %get frame formatting
                    frame = getframe(gcf);
                    %write video frame
                    writeVideo(video,frame);
                    
                end
                
                %close video (otherwise matlab keeps using it and .avi cannot be opened)
                close(video);
            end
            %merge videos
            disp('Merging movies');
            moviename_save=[obj.foldername obj.filename '_trackmovie.avi'];
            video = VideoWriter(moviename_save,'Uncompressed AVI');
            video.FrameRate = 7;
            open(video);
            for qq=1:qqmax
                disp(['Merging movie ' num2str(qq) ' out of ' num2str(qqmax) '.'])
                moviename_load=[obj.foldername obj.filename '_trackmovie_' num2str(qq) '.avi'];
                mov1 = VideoReader(moviename_load);
                while hasFrame(mov1)
                    v=readFrame(mov1);
                    writeVideo(video,v);
                end
            end
            close(video);
            close all;
        end
        function plot_flow_vs_distance_rot(obj)
            %plots the absolute flow velocity versus distance to the
            %droplet in both x and y direction away from the droplet center
            
            radius=120;
            
            if isempty(obj.average_flow_rot)
                error('you need to calculate the average flow profile first.')
            end
            
            %for the y direction in blue
            %midpoint is half length matrix
            xmid=round(length(obj.average_flow_rot.x(1,:))/2);
            yval=obj.average_flow_rot.y(:,xmid);
            speedvaly=sqrt((obj.average_flow_rot.u(:,xmid)).^2+(obj.average_flow_rot.v(:,xmid)).^2);
            
            %for the x direction in red
            %midpoint is half length matrix
            ymid=round(length(obj.average_flow_rot.x(:,1))/2);
            xval=obj.average_flow_rot.y(:,xmid);
            speedvalx=sqrt((obj.average_flow_rot.u(ymid,:)).^2+(obj.average_flow_rot.v(ymid,:)).^2);
            
            figure
            box on
            hold on
            plot(yval*obj.scale/radius,speedvaly*obj.scale*obj.framerate,'bo')
            plot(xval*obj.scale/radius,speedvalx*obj.scale*obj.framerate,'ro')
            xlabel('r (\mum)')
            ylabel('U (\mum/s)')
            set(gca, 'XScale', 'log')
            set(gca, 'YScale', 'log')
            axis([0.8 20 10^-1 10^3])
            hold off
            
        end
        function plot_flow_vs_distance_rotmintrans(obj)
            %plots the absolute flow velocity versus distance to the
            %droplet in both x and y direction away from the droplet center
            
            radius=31;
            
            if isempty(obj.average_flow_rotmintrans)
                error('you need to calculate the average flow profile first.')
            end
            
            %for the y direction in blue
            %midpoint is half length matrix
            xmid=round(length(obj.average_flow_rotmintrans.x(1,:))/2);
            yval=obj.average_flow_rotmintrans.y(:,xmid);
            speedvaly=sqrt((obj.average_flow_rotmintrans.u(:,xmid)).^2+(obj.average_flow_rotmintrans.v(:,xmid)).^2);
            
            %for the x direction in red
            %midpoint is half length matrix
            ymid=round(length(obj.average_flow_rotmintrans.x(:,1))/2);
            xval=obj.average_flow_rotmintrans.y(:,xmid);
            speedvalx=sqrt((obj.average_flow_rotmintrans.u(ymid,:)).^2+(obj.average_flow_rotmintrans.v(ymid,:)).^2);
            
            figure
            box on
            hold on
            plot(abs(yval)*obj.scale/radius,speedvaly*obj.scale*obj.framerate,'bx')
            plot(abs(xval)*obj.scale/radius,speedvalx*obj.scale*obj.framerate,'rx')
            xlabel('r (\mum)')
            ylabel('U (\mum/s)')
            set(gca, 'XScale', 'log')
            set(gca, 'YScale', 'log')
            axis([0.8 10 10^-1 10^3])
            hold off
            
        end
        function obj=output_PIV_mask(obj)
            %npoints indicates the number of points used to draw each mask
            %default is 20;
            
            if nargin<2
                npoints=20;
            end
            
            disp(['Your object contains ' num2str(obj.NOP) ' particles']);
            maskiererx=cell(obj.NOP,obj.NOF);
            maskierery=cell(obj.NOP,obj.NOF);
            for p=1:obj.NOP
                %get user input for radius
                radius(p)=get_numerical_input(['What is the radius you want to use for particle ' num2str(p) '?  ']);
                
                %get trajectory of particle p
                tr_n=NaN(obj.NOF,4);
                tr_n(:,1:2)=obj.tr(find(obj.tr(:,4)==p),1:2);
                tr_n(:,3)=radius(p)*ones(1,obj.NOF);
                tr_n(:,4)=obj.tr(find(obj.tr(:,4)==p),3);
                
                %for each frame, make a mask
                for t=1:obj.NOF
                    x=tr_n(t,1);
                    y=tr_n(t,2);
                    r=tr_n(t,3);
                    mask_tp_x=NaN(npoints,1);
                    mask_tp_y=NaN(npoints,1);
                    for j=1:npoints
                        mask_tp_x(j,1)=x+cos(2*pi*j/npoints)*r;
                        mask_tp_y(j,1)=y+sin(2*pi*j/npoints)*r;
                    end
                    maskiererx{p,t}=mask_tp_x;
                    maskierery{p,t}=mask_tp_y;
                end
            end
            obj.contourfile={maskiererx,maskierery};
            %save the produced mask in the current folder
            save('contour.mat', 'maskiererx', 'maskierery');
            
        end%end function output PIV
    end %End methods
end %End object