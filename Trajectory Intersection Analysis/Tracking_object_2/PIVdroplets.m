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
    properties
        Flows %1xNOF cell structure that contains a profile_class object in each cell
        average_flow %1x1 structure of type profileclass containing the average flow profile
        average_flow_trans %1x1 structure of type profileclass containing the average flow profile after only translation
        average_flow_rot %1x1 structure of type profileclass containing the average flow profile after translation and rotation
        average_flow_rotmintrans %1x1 structure of type profileclass containing the average rotated profile with the average translated profile substracted
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
            for n=2:obj.NOF
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
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            %Make sure flows are present;
            if isempty(obj.Flows)
                obj=add_PIVdata(obj);
            end
            
            % Get particle trace
            trace=obj.tr;
            
            %Check if selection option is used
            if selection==0
                nframes=length(trace(:,1))-1;
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
                    X(k,:)=trace(k,1);
                    Y(k,:)=trace(k,2);
                end
            end
            
            % Create list of translated and rotated profiles
            profile_list=[];
            for k=1:nframes
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
                profile_list=[profile_list framek];
            end
            
            % Average and plot average
            averaged_profile=average_profiles_fast(profile_list);
            %plot_profile(averaged_profile)
            
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
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            %Make sure flows are present;
            if isempty(obj.Flows)
                obj=add_PIVdata(obj);
            end
            
            % Get particle trace
            trace=obj.tr;
            
            %Check if selection option is used
            if selection==0
                nframes=length(trace(:,1))-1;
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
                    X(k,:)=trace(k,1);
                    Y(k,:)=trace(k,2);
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
            averageflow=translate_profile(averageflow,imsize(2),imsize(1));
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
                       
            %plot object
            figure
            hold on
            if imcheck~=0
                colormap('gray'), imshow(image);
            end
            box on
            axis([0 max(max(x)) 0 max(max(y))]);
            pbaspect([max(max(x))/max(max(y)) 1 1]);
            quiver(x,y,u,v,1,'b');
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
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
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
            nframes=obj.NOF;
            
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
            trace=obj.tr;
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
                    Ik=imread(fname,k);
                    
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
    end %End methods
end %End object