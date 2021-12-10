classdef particles
    %Particles class
    %Takes as input an 8-bit grayscale tiffstack of colloidal particles in 2D.
    %Use to find and track particles in image and subsequently calculate
    %various properties of the system
    %All calculated values are saved as an object property
    %
    %Available functions:
    % - obj=particles()                 Initialization, gets filename
    % - obj=add_pos(obj,pos)            manually add particle position matrix
    % - obj=add_trace(obj,trace)        manually add particle trajectory matrix
    % - obj=add_image_stack(obj,A)      manually add image stack
    % - output=get_parameters(obj)      returns parameters used for finding and tracking particles
    % - output=get_name(obj)            returns the location of the imagestack connected to this object
    % - obj=read_parameters(obj)        reads in parameters for particle finding and Ftracking. Uses format that comes out of get_parameters.
    % - obj=find_pos(obj)               finds particle positions in each frame. Interactive.
    % - check_pixel_bias(obj)           plots a graph that allows you to check for pixel bias
    % - output=get_positions(obj)       returns a list of all particle positions
    % - obj=change_scale(obj,scale)     saves a new scale to object in micrometer per pixel
    % - obj=change_framerate(obj)       saves a new framerate to object in frames per second
    % - obj=track_obj(obj)              Uses Crocker and Grier code to connect particle positions into trajectories
    % - output=get_trace(obj)           returns a list of all paritcle trajectories
    % - plot_tracks(obj,tbegin,tend,ID) plots particle trajectories
    % - displacement_histogram(obj)     plots displacement histograms
    % - obj=drift_correct(obj)          corrects for constant drift in the sample
    % - obj=make_dist_matrix(obj)       creates a table of all average distances between particles
    % - obj=make_av_dist_matrix(obj)    creates a table of the average distances between all particles
    % - plot_trace_on_movie(obj)        makes movie of particle motion with arrow of direction of motion on top
    % - output=find_abs_displacement(obj)       returns the particles that have a larger displacement than 3 times over the average
    % - obj=complete_trace(obj)         fills holes in trajectories with NaN
    % - obj=automatically_improve_trace(obj)    uses find_abs_displacement to find likely mistracks, removes them and tracks again
    % - plot_tracked_image(obj,frame)   plots image with vector of direction of motion on top
    % - obj=remove_point_trace(obj)     manually removes one frame from the trajectories
    % - obj=merge_tracks(obj,ID1,ID2)   Merges trajectories with label 1 and label 2 if they don't overlap
    % - obj=calc_MSD(obj)               calculates the mean squared displacement
    % - plot_MSD(obj)                   plots the mean squared displacement vs time for each particle and the average
    % - output=find_stagnant_particles(obj)     returns ID's of particles that move less than 3 times below the average
    % - obj=remove_stagnant_particles(obj)      uses find_stagnant_particles to remove particles that are presumed stuck
    % - obj=find_reduced_trajectory(obj)        calculates a 2D histogram of direction and magnitude of steps
    % - obj=plot_reduced_tracks(obj)            plots the histogram calculated by find_reduced_trajectory
    % - obj=plot_average_reduced_track(obj)     see two functions above: plotst the average histogram
    % - obj=plot_average_reduced_track_subplot(obj,dt)  plots the histogram for different time intervals dt and places in subplot
    % - output=find_reduced_trajectory_temp(obj,dt)     returns the reduced trajectory
    % - obj=correlate_motion(obj)       uses reduced tracks to calculte Dtrans, Drot and v0(persistent speed)
    % - obj=crop_in_time(obj,tmin,tmax)         removes all information outside of window tmin tmax
    % - obj=calc_speed_vs_tiume(obj,particlelist)       for all particles in the list, calculates speed vs time
    % - obj=calc_speed_vs_time_neighbours(obj, ID)      correlates speed of particle with ID with that of neighbours
    % - obj=find_neighbours(obj,cutoff)         finds neighbours of all particles using distance matrix. All within cutoff are neighbours
    % - obj=find_neighbours_time(obj,cutoff)    finds number of neighbours in each frame
    % - correlate_speed_numneighs(obj,cutoff,particlelist)  correlates particle speed with number of neighbours
    % - output=corr_xy_motion(obj)      finds average diffusion tensor for all particles in the object
    %
    %
    %Requires the following external files to operate
    % - bpass.m (by Crocker and Grier)
    % - ctrd.m (by Crocker and Grier)
    % - pkfnd.m (by Crocker and Grier)
    % - track.m (by Crocker and Grier)
    % - correlate_functions.m
    % - distance_filter.m
    % - findchain_class.m
    % - get_numerical_input.m
    % - get_textual_input.m
    % - initial_particle_peakfind_class.m
    % - peakfindoptimizer_class.m
    % - secondary_particle_peakfind_class.m
    %
    %Make sure all of the above functions are located in the same folder as
    %this code.
    %
    %Created on 19-07-17 by Pepijn Moerman
    %Dates modified: 06-10-2017
    %                14-12-2017
    %                26-02-2018
    %
    % Note: to save produced images as vector format use:
    % fig1=figure(1);
    % fig1.Renderer='Painters';
    % alternatively use export_fig('name','painters')
    
    properties
        foldername      %name of folder in which tiffstack is located. Saves data to this folder also
        filename        %identification name of this object
        scale           %scale of the image in um/pixel
        framerate       %framerate of the movie
        pos             %positions of particles
        tr              %n x 4 matrix containing all particle positions in time and a label
        motion          %2 x 1 matrix containig velocity and diffusion coefficient of particles
        gr              %radial distribution function for each frame
        MSD             %Mean Squared Displacement of all particles individually
        VCF             %velocity correlation function of all individual particles
        VCFav           %average velocity correlation function for all particles in the sample
        speeds          %matrix of width NOP and length NOF containing instantaneous speed at each time for each particle
        disps           %cell containing dx displacement and dy displacement for objects;
        neighbours      %list containing which particles are neighbours
        num_neighs      %average number of neighours
        numneighs_time  %number of neighbours in each frame
        cutoff          %cutoff used for calculating num_neighs
        cutoff_time     %cutoff used for calculating numneighs_time
        tr_red          %Reduced trajectory
        im              %Image of particles in frame 1
        dist            %distance matrix
        dist_av         %average distance matrix
        drift_corr_hist %history of drift corrections
        ydriftval       %average drift in y direction: for gradient experimetns
        posparam        %parameters used for finding particles
        trparam         %parameters used for tracking particles
        NOP             %number of particles
        NOF             %number of frames
        
    end
    methods
        function obj=particles(fname,parameters)
            %Takes n x 3 matrix pos as input and saves in object
            %Pos should be organized as [x y f]
            %Here x and y are coordintes and f is framenumber
            %Parameters needs to be a 4x2 cell obtained using the
            %get_parameters function from this object;
            
            if nargin==2
                if isempty(fname)
                    disp('Select folder that contains the tiffstack(s)')
                    [filen, foldername_temp] = uigetfile('*.tif','Select tif stack');
                    obj.scale=1;
                    obj.framerate=1;
                else
                    startIndex=max(regexp(fname,'\'));
                    foldername_temp=fname(1:startIndex);
                    filen=fname(startIndex+1:end);
                end
                try obj=read_parameters(obj,parameters);
                catch
                    disp('One or more parameters were not found.')
                end
            elseif nargin==1
                startIndex=max(regexp(fname,'\'));
                foldername_temp=fname(1:startIndex);
                filen=fname(startIndex+1:end);
                obj.scale=1;
                obj.framerate=1;
            else
                disp('Select folder that contains the tiffstack(s)')
                [filen, foldername_temp] = uigetfile('*.tif','Select tif stack');
                obj.scale=1;
                obj.framerate=1;
            end
            
            try A = imread([foldername_temp filen],1);
            catch
                disp(['No tifstack exists with name: ' foldername_temp filen '.'])
                disp('The imagestack is now "NaN". Be aware that some functions might not work')
                A=NaN;
                foldername_temp=NaN;
                filen=NaN;
            end
            obj.im=A;
            obj.foldername=foldername_temp;
            if isnan(filen)
                obj.filename=NaN;
            else
                obj.filename=regexprep(filen,'.tif','');
            end
            
        end
        function obj=find_pos_hough(obj,initialize)
            %Creates obj.pos by tracking particles using hough transform
            %Function is interactive and requires human input.
            %Initialize indicates whether or not existing parameters are
            %used for particle_finding or a new set is generated
            %interactively. Use existing set is (0), generate new set
            %interactively is (1). Default is (1)
            
            %set off warning
            id='images:imfindcircles:warnForLargeRadiusRange';
            warning('off',id)
            
            if nargin<2
                initialize=1;
            end
            if isempty(obj.posparam)
                initialize=1;
            end
            fname=[obj.foldername obj.filename '.tif'];
            if initialize==1
                output=dropletfinder_class(fname);
            else
                output=dropletfinder_class(fname,obj.posparam);
            end
            obj.pos=output{1};
            obj.posparam=output{2};
            
        end
        function obj=add_pos(obj,pos_temp)
            %Takes n x 3 matrix pos_temp as input and saves in object
            %Pos should be organized as [x y f]
            %Here x and y are coordintes and f is framenumber
            obj.pos=pos_temp(:,1:3);
        end
        function obj=add_trace(obj,tr_temp)
            %Takes n x 4 matrix tr_temp as input and saves in object
            %Pos should be organized as [x y f]
            %Here x and y are coordintes and f is framenumber
            obj.tr=tr_temp(:,1:4);
            obj.NOP=max(tr_temp(:,4));
            obj.NOF=max(tr_temp(:,3));
            obj=complete_trace(obj);
            
        end
        function obj=add_image_stack(obj,A_temp)
            %Manually replaces the object image_stack obj.im with another
            %imagestack A_temp
            obj.im=A_temp;
        end
        function parameters=get_parameters(obj)
            %outputs a cell containing the parameters used for tracking in
            %this object
            parameters=cell(4,2);
            parameters{1,1}='framerate';
            parameters{1,2}=obj.framerate;
            parameters{2,1}='scale';
            parameters{2,2}=obj.scale;
            parameters{3,1}='position finding parameters';
            parameters{3,2}=obj.posparam;
            parameters{4,1}='tracking parameters';
            parameters{4,2}=obj.trparam;
        end
        function fname=get_name(obj)
            %Outputs the name of the folder and file where the data are located.
            %Use to create new object with same file
            fname=[obj.foldername obj.filename '.tif'];
        end
        function obj=read_parameters(obj,parameters)
            %parameters needs to be a 4x2 cell created using the
            %get_parameters function of this object
            obj.framerate=parameters{1,2};
            obj.scale=parameters{2,2};
            obj.posparam=parameters{3,2};
            obj.trparam=parameters{4,2};
        end
        function obj=find_pos(obj,initialize)
            %Creates obj.pos by tracking particles
            %Function is interactive and requires human input
            %initialize indicates whether or not existing parameters are
            %used for particle_finding or a new set is generated
            %interactively. Use existing set is (0), generate new set
            %interactively is (1). Default is (1)
            if nargin<2
                initialize=1;
            end
            if isempty(obj.posparam)
                initialize=1;
            end
            MSGID3='MATLAB:imagesci:tifftagsread:expectedTagDataFormatMultiple';
            warning('off', MSGID3)
            fname=[obj.foldername obj.filename '.tif'];
            if initialize==1
                output=particlefinder_class(fname);
            else
                output=particlefinder_class(fname,obj.posparam);
            end
            obj.pos=output{1};
            obj.posparam=output{2};
            
            %check pixel bias
            check_pixel_bias(obj)
            
        end
        function obj=find_pos_selection(obj,n,initialize)
            %Creates obj.pos by tracking particles
            %Function is interactive and requires human input
            %initialize indicates whether or not existing parameters are
            %used for particle_finding or a new set is generated
            %interactively. Use existing set is (0), generate new set
            %interactively is (1). Default is (1)
            if nargin<3
                initialize=1;
            end
            if nargin<2
                error('You need to specify the number of images you want to analyze. If you want to analyze all, use find_pos(obj).')
            end
            if isempty(obj.posparam)
                initialize=1;
            end
            MSGID3='MATLAB:imagesci:tifftagsread:expectedTagDataFormatMultiple';
            warning('off', MSGID3)
            fname=[obj.foldername obj.filename '.tif'];
            if initialize==1
                output=particlefinder_class_selection(fname,n);
            else
                output=particlefinder_class_selection(fname,n,obj.posparam);
            end
            obj.pos=output{1};
            obj.posparam=output{2};
            
            %check pixel bias
            check_pixel_bias(obj)
            
        end
        function obj=find_pos_confined(obj,n,initialize)
            %Creates obj.pos by tracking particles
            %Function is interactive and requires human input
            %initialize indicates whether or not existing parameters are
            %used for particle_finding or a new set is generated
            %interactively.
            %Tailored for particles close to boundaries
            %n is the number of gray pixels you want to add to each
            %boundary
            if nargin<2
                n=20;
            end
            if nargin<3
                initialize=1;
            end
            if isempty(obj.trparam)
                initialize=1;
            end
            MSGID3='MATLAB:imagesci:tifftagsread:expectedTagDataFormatMultiple';
            warning('off', MSGID3)
            fname=[obj.foldername obj.filename '.tif'];
            if initialize==1
                output=particlefinder_class2(fname,n);
            else
                output=particlefinder_class2(fname,n,obj.posparam);
            end
            obj.pos=output{1};
            obj.pos(:,1:2)=obj.pos(:,1:2)-n;
            obj.posparam=output{2};
            
            %check pixel bias
            check_pixel_bias(obj)
            
        end
        function check_pixel_bias(obj)
            %checks for pixel bias by plotting residue of x and y
            %coordinates
            
            if isempty(obj.pos)
                disp('Why are you checking for pixel bias? You have not even started finding the particle positions yet?')
                obj=find_pos(obj);
            end
            posi=obj.pos;
            x_res=mod(posi(:,1),1);
            y_res=mod(posi(:,2),1);
            figure
            hold on
            box on
            histogram(x_res);
            histogram(y_res);
            xlabel('Residual coordinate (pixel)')
            ylabel('P(residual coordinate)')
        end
        
        function pos=get_positions(obj)
            %Gets the positions from object
            pos=obj.pos;
        end
        function obj=change_scale(obj,scale_temp)
            %Manually change the scale
            obj.scale=scale_temp;
        end
        function obj=change_framerate(obj,framerate_temp)
            %Manually change the framerate
            obj.framerate=framerate_temp;
        end
        function obj=track_obj(obj,initialize)
            %Takes obj and tracks the particle positions
            %outputs a n x 4 matrix organized as [x y f ID]
            %Initialize indicates whether or not to use existing parameters
            %for tracking or define new ones. Use existing ones is (0) and
            %define new ones is (1). Default is (1).
            if nargin<2
                initialize=1;
            end
            if isempty(obj.trparam)
                initialize=1;
            end
            
            if isempty(obj.pos)
                obj=find_pos(obj);
            end
            
            if initialize==1
                % Set parameters for trackingsoftware
                maxdisp=get_numerical_input('What is the max displacement? ');
                param.mem=get_numerical_input('How long can a particle go missing? ');
                param.dim=2; % dimensionality of data
                param.good=get_numerical_input('What is the minimal number of frames for a track to be accepted? ');
                param.quiet=1; % 0 = text, 1 = no text
            else
                maxdisp=obj.trparam{2};
                param=obj.trparam{1};
            end
            
            moveon=1;
            try trace=track(obj.pos,maxdisp,param);
            catch
                disp('Tracking failed. Try using different parameters');
                moveon=0;
            end
            if moveon==1
                obj.NOP=max(trace(:,4));
                obj.NOF=max(trace(:,3));
                obj.tr=trace;
                obj.trparam={};
                obj.trparam{1}=param;
                obj.trparam{2}=maxdisp;
                obj=complete_trace(obj);
            end
        end
        function obj=track_obj_ydrift(obj,factor,initialize)
            %Connects particle positions of a movie with large y-drift by
            %first decreasing the y coordinates by a factor 'factor', then
            %tracking and finally multiplying the y coordinates back to the
            %original.
            
            obj.ydriftval=factor;
            if nargin<3
                initialize=1;
            end
            if isempty(obj.trparam)
                initialize=1;
            end
            
            if isempty(obj.pos)
                obj=find_pos(obj);
            end
            
            pos_temp=[obj.pos(:,1)/factor obj.pos(:,2) obj.pos(:,3)];
            
            if initialize==1
                % Set parameters for trackingsoftware
                maxdisp=get_numerical_input('What is the max displacement? ');
                param.mem=get_numerical_input('How long can a particle go missing? ');
                param.dim=2; % dimensionality of data
                param.good=get_numerical_input('What is the minimal number of frames for a track to be accepted? ');
                param.quiet=1; % 0 = text, 1 = no text
            else
                maxdisp=obj.trparam{2};
                param=obj.trparam{1};
            end
            
            moveon=1;
            try trace=track(pos_temp,maxdisp,param);
            catch
                disp('Tracking failed. Try using different parameters');
                moveon=0;
            end
            if moveon==1
                trace=[trace(:,1)*factor trace(:,2) trace(:,3:4)];
                obj.NOP=max(trace(:,4));
                obj.NOF=max(trace(:,3));
                obj.tr=trace;
                obj.trparam={};
                obj.trparam{1}=param;
                obj.trparam{2}=maxdisp;
                obj=complete_trace(obj);
            end
            
        end
        function tr=get_trace(obj)
            %Gets the traces from object
            tr=obj.tr;
        end
        function plot_tracks(obj,begint,endt,ID)
            %Plots the particle trajectories in time over the first image
            %of the move
            %If you want to only look at trajectory ID, but all frames,
            %put 0 for begint and endt.
            
            %             if isnan(obj.im)
            %                 error('There is no imagestack loaded in the object, so no trajectory can be plotted')
            %             end
            
            %if trace does not exist, trace particles first
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            
            if nargin<3
                endt=max(obj.tr(:,3));
                begint=min(obj.tr(:,3));
                ID=0;
            elseif nargin<2
                begint=min(obj.tr(:,3));
                ID=0;
            elseif nargin==3
                ID=0;
                if begint==0 && endt==0
                    begint=min(obj.tr(:,3));
                end
            elseif nargin==4 && begint==0 && endt==0
                begint=min(obj.tr(:,3));
                endt=obj.NOF;
            end
            
            if any(ID>obj.NOP)
                disp('The given ID does not exist. Please try a different value')
            else
                
                traj=get_trace(obj);
                if ~isnan(obj.im)
                    A=imread([obj.foldername obj.filename '.tif'],begint);
                    %picture=obj.im;
                    picture=A;
                    %Assert image is 8-bit
                    if length(size(picture))>2
                        %Convert to grayscale
                        picture = rgb2gray(picture);
                    end
                    picturesize=size(obj.im);
                    size1=picturesize(1);
                    size2=picturesize(2);
                end
                
                
                cmap=jet(obj.NOP);
                figure
                hold on
                box on
                if ~isnan(obj.im)
                    colormap('gray'),imagesc(picture);
                    %axis([0 length(picture(1,:)) 0 length(picture(:,1))]);
                end
                if any(ID==0)
                    for j=1:max(traj(:,4))
                        trj=traj(find(traj(:,4)==j),1:3);
                        plot(trj(begint:endt,1),trj(begint:endt,2),'Color',cmap(j,:),'DisplayName',['Particle ' num2str(j)],'LineWidth',1);
                    end
                else
                    for j=1:max(traj(:,4))
                        if any(ID==j)
                            trj=traj(find(traj(:,4)==j),1:3);
                            plot(trj(begint:endt,1),trj(begint:endt,2),'Color',cmap(j,:),'DisplayName',['Particle ' num2str(j)],'LineWidth',1);
                        end
                    end
                end
                xlabel('X (pixels)','FontSize',18)
                ylabel('Y (pixels)','FontSize',18)
                lgd.Color=[0.5 0.5 0.5];
                lgd.TextColor=[1 1 1];
                lgd=legend('show','Location','northeastoutside');
                pbaspect([size2 size1 1])
                axis([0 size2 0 size1])
                hold off
            end
        end
        function plot_tracks_selection(obj,fraction)
            %Plots the particle trajectories in time over the first image
            %of the move
            %If you want to only look at trajectory ID, but all frames,
            %put 0 for begint and endt.
            
            %             if isnan(obj.im)
            %                 error('There is no imagestack loaded in the object, so no trajectory can be plotted')
            %             end
            
            %if trace does not exist, trace particles first
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            
            if nargin<2
                fraction=1;
            elseif fraction>1
                fraction=1;
            elseif fraction<0
                fraction=0;
            end
            traj=obj.tr;
            
            figure
            hold on
            box on
            if fraction>0
                step=floor(1/fraction);
            else
                step=floor(obj.NOP/2);
            end
            if step<1
                step=1;
            end
            step
            for j=1:step:obj.NOP
                trj=traj(find(traj(:,4)==j),1:3);
                plot(trj(:,1),trj(:,2),'Color','r','LineWidth',1);
            end
            xlabel('X (pixels)','FontSize',18)
            ylabel('Y (pixels)','FontSize',18)
            pbaspect([1 1 1])
            axis([0 512 0 512])
            hold off
            
        end
        function plot_tracks_zero(obj,lengthdata,IDlist)
            %Draws trajectories like Howse paper
            %First translates such that t=0 starts at x,y=0
            %Plots the particle trajectories in time
            
            
            %if trace does not exist, trace particles first
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            
            %make sure enough input is provided
            if nargin<3
                error('Please provide all input parameters.')
            end
            
            %make sure lengthdata and IDlist are acceptable input
            if any(IDlist>obj.NOP)
                disp('The given ID does not exist. Please try a different value')
            end
            if lengthdata>obj.NOF
                error('Not all trajectories are long enough. Please provide a shorter length.')
            end
            
            
            traj=get_trace(obj);
            
            cmap=[0 0 0; 0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0];
            figure
            hold on
            box on
            axis([-700 700 -700 700])
            pbaspect([1 1 1])
            count=0;
            for p=1:obj.NOP
                if any(IDlist==p)
                    count=count+1;
                    trp=traj(find(traj(:,4)==p),1:2);
                    trp=trp(~isnan(trp(:,1)),1:2);
                    trp(:,1)=trp(:,1)-trp(1,1);
                    trp(:,2)=trp(:,2)-trp(1,2);
                    try plot(trp(1:lengthdata,1),trp(1:lengthdata,2),'Color',cmap(count,:),'LineWidth',2);
                    catch
                        error(['Trajectory of particle ' num2str(p) ' was not long enough. Use a shorter length or different IDs'])
                    end
                end
            end
            hold off
            
        end
        function plot_tracks_pretty(obj,begint,endt,ID)
            %Plots the particle trajectories in time over the first image
            %of the move
            %If you want to only look at trajectory ID, but all frames,
            %put 0 for begint and endt.
            
            %             if isnan(obj.im)
            %                 error('There is no imagestack loaded in the object, so no trajectory can be plotted')
            %             end
            
            %if trace does not exist, trace particles first
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            
            if nargin<3
                endt=obj.NOF;
                begint=min(obj.tr(:,3));
                ID=1;
            elseif nargin<2
                begint=min(obj.tr(:,3));
                ID=1;
            elseif nargin==3
                ID=1;
                if begint==0 && endt==0
                    begint=min(obj.tr(:,3));
                end
            elseif nargin==4 && begint==0 && endt==0
                begint=min(obj.tr(:,3));
                endt=obj.NOF;
            end
            
            if any(ID>obj.NOP)
                disp('The given ID does not exist. Please try a different value')
            else
                
                traj=get_trace(obj);
                if ~isnan(obj.im)
                    picture=obj.im;
                    %Assert image is 8-bit
                    if length(size(picture))>2
                        %Convert to grayscale
                        picture = rgb2gray(picture);
                    end
                end
                
                %method 1
                %n=(endt-begint+1);
                
                %method 2
                %get n for the colormap
                n=(endt-begint);
                
                %get trajectory of particle ID
                trj=traj(find(traj(:,4)==ID),1:3);
                %make sure trajectory contains no NAN's, but goes from
                %begint to endt otherwise
                begint_temp=min(trj(~isnan(trj(:,1)),3));
                endt_temp=max(trj(~isnan(trj(:,1)),3));
                endt=min([endt endt_temp]);
                begint=max([begint begint_temp]);
                
                %get n for selection colormap in case there were NaN's
                %important to keep color code constant
                
                %method 1
                %                 nnew=endt-begint+1;
                
                %method 2
                nnew=endt-begint;
                
                figure
                hold on
                box on
                %                 if ~isnan(obj.im)
                %                     colormap('gray'),imagesc(picture);
                %                     axis([0 length(picture(1,:)) 0 length(picture(:,1))]);
                %                 end
                %plot the trajectory
                x=trj(begint:endt,1).';
                y=trj(begint:endt,2).';
                axis equal
                size(x);
                size(y);
                %method 1
                %                 p=plot(x,y,'r','LineWidth',3);
                %                 % modified jet-colormap
                %                 cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
                %                 %take only the selection of cd you want to plot
                %                 cd = cd(:,1:nnew);
                %                 size(cd)
                %                 drawnow
                %                 set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
                
                %method 2
                cmap=jet(n);
                for i=1:nnew
                    plot(x(i:i+1),y(i:i+1),'Color',cmap(i,:),'LineWidth',2)
                end
                axis([0 1200 000 1200])
                hold off
                
            end
            
        end
        function plot_xy_in_time(obj)
            for p=1:obj.NOP
                tr1=obj.tr(find(obj.tr(:,4)==p),1:3);
                figure
                hold on
                box on
                plot(tr1(:,3)/16,(tr1(:,2)-tr1(1,2))*1.44,'r.-','MarkerSize',10)
                plot(tr1(:,3)/16,(tr1(:,1)-tr1(1,1))*1.44,'b.-','MarkerSize',10)
                hold off
            end
        end
        function plot_tracks_drift_corr(obj)
            %Plot tracks, but not on movie but in space. Makes sence only
            %after a drift correction. Useful because some tracks fall off
            %the image after a drift correction.
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            
            
            traj=get_trace(obj);
            cmap=jet(max(traj(:,4)));
            figure
            hold on
            box on
            for j=1:max(traj(:,4))
                trj=traj(find(traj(:,4)==j),1:3);
                plot(trj(:,1),trj(:,2),'Color',cmap(j,:),'DisplayName',['Particle ' num2str(j)]);
            end
            axis([0.9*min(traj(:,1)) 1.1*max(traj(:,1)) 0.9*min(traj(:,2)) 1.1*max(traj(:,2))])
            hold off
        end
        function plot_xy_speed(obj,ID,dt)
            %plots the speed in x and y direction seperately
            %works for only one particle at a time. If no ID is provided,
            %uses the first particle
            
            if nargin<2
                ID=1;
            end
            if nargin<3
                dt=1;
            end
            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            trj=obj.tr(find(obj.tr(:,4)==ID),1:3);
            x=trj(:,1);
            x=smooth(x);
            y=trj(:,2);
            y=smooth(y);
            t=trj(:,3);
            dx=x(1:dt:end-dt,1)-x(1+dt:dt:end);
            dy=y(1:dt:end-dt,1)-y(1+dt:dt:end);
            
            figure
            hold on
            box on
            %plot(t(1:dt:end-1)/obj.framerate,dx*obj.scale*obj.framerate,'ro');
            %plot(t(1:dt:end-dt)/obj.framerate,-dy*obj.scale*obj.framerate,'ko');
            plot((t(57:125)-57)/obj.framerate,-dy(57:125,1)*obj.scale*obj.framerate,'ko');
            axis([0 15 -700 1000])
            xticks([0 5 10 15])
            yticks([-500 0 500 1000])
            hold off
            
            
        end
        function displacement_histogram(obj,nbins,show_individuals)
            %Plots histograms of instantaneous displacement in x and y
            %Uses the tracked coordinates and plots the average
            %displacement histogram
            %show individuals determines whether or not histograms for
            %individual cases are shown. Pick (1) for yes and (0) for no.
            %Default is (0);
            if nargin<2
                nbins=100;
            end
            if nargin<3
                show_individuals=0;
            end
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            trace=get_trace(obj);
            nframes=obj.NOF;
            nparts=obj.NOP;
            cmap=jet(nparts);
            dx_tot=[];
            dy_tot=[];
            if show_individuals==1
                figure
                title('Displacements dx')
                hold on
                box on
                lgd=legend('show');
            end
            for p=1:nparts
                trj=trace((p-1)*nframes+1:p*nframes,1:2);
                x1=trj(1:length(trj(:,1))-1,1);
                x2=trj(2:length(trj(:,1)),1);
                dx=x2-x1;
                dx_tot=[dx_tot ; dx];
                if show_individuals==1
                    histogram(dx,'FaceColor',cmap(p,:),'DisplayName',['Particle ' num2str(p)])
                end
            end
            if show_individuals==1
                xlabel('Displacement (pixels)','FontSize',18)
                ylabel('Number','FontSize',18)
                lgd.Color=[1 1 1];
                lgd.TextColor=[0 0 0];
                hold off
                
                figure
                title('Displacements dy')
                hold on
                box on
                lgd=legend('show');
            end
            for p=1:nparts
                trj=trace((p-1)*nframes+1:p*nframes,1:2);
                y1=trj(1:length(trj(:,1))-1,2);
                y2=trj(2:length(trj(:,1)),2);
                dy=y2-y1;
                dy_tot=[dy_tot ; dy];
                if show_individuals==1
                    histogram(dy*obj.scale,'FaceColor',cmap(p,:),'DisplayName',['Particle ' num2str(p)])
                end
            end
            if show_individuals==1
                xlabel('Displacement (\mum)','FontSize',16)
                ylabel('Number','FontSize',16)
                lgd.Color=[1 1 1];
                lgd.TextColor=[0 0 0];
                hold off
            end
            
            figure
            hold on
            box on
            title('Average displacement in x and y')
            %histogram(dx_tot*obj.scale)
            %histogram(dy_tot*obj.scale)
            h=histfit(dx_tot*obj.scale,nbins);
            h(1).FaceColor='b';
            h(1).FaceAlpha=0.5;
            h(2).Color='b';
            %h(2).LineSpec='- -';
            g=histfit(dy_tot*obj.scale,nbins);
            g(1).FaceColor='r';
            g(1).FaceAlpha=0.5;
            g(2).Color='r';
            legend('','x','','y')
            xlabel('Displacement (\mum)','FontSize',16)
            ylabel('PDF','FontSize',16)
            hold off
            
            figure
            hold on
            box on
            title('Average logarithmic displacement in x and y')
            h=histfit(dx_tot*obj.scale,nbins);
            h(1).FaceColor='b';
            h(1).FaceAlpha=0.5;
            h(2).Color='b';
            %h(2).LineSpec='- -';
            g=histfit(dy_tot*obj.scale,nbins);
            g(1).FaceColor='r';
            g(1).FaceAlpha=0.5;
            g(2).Color='r';
            %g(2).LineSpec='- -';
            set(gca,'YScale','log')
            legend('','x','','y')
            xlabel('Displacement (\mum)','FontSize',16)
            ylabel('PDF','FontSize',16)
            hold off
            
        end
        function bin_disp_hist(obj,nbins,driftval)
            %Splits the image in nbins number of bins in the y dimension
            %Default number of bins is 7
            %Driftval is the drift correction. Use number of pixels a
            %particle moves on average in the x direction for this
            
            if nargin<2
                nbins=7;
            end
            if nargin<3
                driftval=5;
            end
            
            if isempty(obj.pos)
                'wtf?'
                obj=find_pos(obj);
            end
            
            if ~isempty(obj.im)
                channel_height=length(obj.im(:,1));
            else
                error('An image must be provided for this function to work.')
            end
            
            value=NaN(nbins,2);
            
            for b=1:nbins
                binlim_low=channel_height/nbins*(b-1);
                binlim_up=channel_height/nbins*b;
                pos_b=obj.pos(find(obj.pos(:,2) <= binlim_up & obj.pos(:,2) >= binlim_low),:);
                if b==1
                    obj_b1=particles(get_name(obj));
                    obj_b1=add_pos(obj_b1,pos_b);
                    obj_b1=track_obj_ydrift(obj_b1,driftval);
                    obj_b=obj_b1;
                else
                    obj_b=particles(get_name(obj),get_parameters(obj_b1));
                    obj_b=add_pos(obj_b,pos_b);
                    obj_b=track_obj_ydrift(obj_b,driftval,0);
                end
                obj_b=remove_particles_flow(obj_b,driftval);
                try value(b,2)=disp_hist_xy_sep(obj_b);
                catch
                    value(b,2)=NaN;
                end
                value(b,1)=(b-0.5)*channel_height/nbins;
            end
            
            figure
            hold on
            box on
            plot(value(:,1),value(:,2),'ko')
            hold off
            
        end
        function obj=disp_hist_xy_sep(obj,nbins,shiftval)
            %Display lin-lin displacmenet histograms of separate x and y
            %displacements
            
            if isempty(obj.tr)
                error('No particles in object');
            end
            
            if nargin<2
                nbins=10;
            end
            
            if nargin<3
                shiftval=0;
            end
            
            trace=get_trace(obj);
            nframes=obj.NOF;
            nparts=obj.NOP;
            dx_tot=[];
            dy_tot=[];
            
            for p=1:nparts
                trj=trace((p-1)*nframes+1:p*nframes,1:2);
                y1=trj(1:length(trj(:,1))-1,2);
                y2=trj(2:length(trj(:,1)),2);
                dy=y2-y1;
                yav=(trj(2:end,2)+trj(1:end-1,2))/2;
                %dy_tot=[dy_tot ; dy yav];
                dy_tot=[dy_tot ; dy];
            end
            dy_mean=nanmean(dy_tot(:,1));
            dy_tot=dy_tot+shiftval;
            dy_mean1=nanmean(dy_tot(1:end/3,1));
            dy_mean2=nanmean(dy_tot(end/3:2*end/3,1));
            dy_mean3=nanmean(dy_tot(2*end/3:end,1));
            %dy_std=nanstd(dy_tot(:,1));
            dy_std=nanstd([dy_mean1 dy_mean2 dy_mean3]);
            disp(['The average displacement in y is ' num2str(dy_mean*obj.scale*obj.framerate) '.'])
            disp(['The error bar in x is ' num2str(dy_std*obj.scale*obj.framerate) '.'])
            
            for p=1:nparts
                trj=trace((p-1)*nframes+1:p*nframes,1:2);
                x1=trj(1:length(trj(:,1))-1,1);
                x2=trj(2:length(trj(:,1)),1);
                dx=x2-x1;
                dx_tot=[dx_tot ; dx];
            end
            dx_mean=nanmean(dx_tot);
            dx_std=nanstd(dx_tot);
            disp(['The average displacement in x is ' num2str(dx_mean*obj.scale*obj.framerate) '.'])
            disp(['The error bar in x is ' num2str(dx_std*obj.scale*obj.framerate) '.'])
            
                        figure
                        hold on
                        box on
                        title('Average displacement in x')
                        h=histfit(dx_tot*obj.scale*obj.framerate,nbins);
                        h(1).FaceColor='b';
                        h(1).FaceAlpha=0.5;
                        h(2).Color='b';
                        xlabel('Displacement (\mum)','FontSize',16)
                        ylabel('PDF','FontSize',16)
                        hold off
            
            
            figure
            hold on
            box on
            title('Average displacement in y')
            g=histfit(dy_tot(:,1)*obj.scale*obj.framerate,nbins);
            g(1).FaceColor='b';
            g(1).FaceAlpha=1;
            g(1).EdgeColor='k';
            g(1).LineWidth=1;
            g(2).Color='r';
            %xlim([-1 1])
            %xticks([]);
            %yticks([]);
            xlabel('Displacement (\mum)','FontSize',16)
            ylabel('PDF','FontSize',16)
            hold off
            
            obj.disps={dx_tot,dx_std,dy_tot,dy_std};
        end
        function obj=analyze_grad_disp(obj,obj_ref,lowerthreshold,upperthreshold)
            
            obj=read_parameters(obj,get_parameters(obj_ref));
            driftval=27;
            obj.ydriftval=driftval;
            if nargin<3
                lowerthreshold=-45;
            end
            if nargin<4
                upperthreshold=-10;
            end
            if isempty(obj.pos)
                obj=find_pos_confined(obj,20,0);
            end
            %if isempty(obj.tr)
                obj=position_distance_selection(obj,20);
                obj=track_obj_ydrift(obj,driftval,0);
                
            %end
            if ~isempty(obj.tr)
                obj=remove_particles_flow(obj,lowerthreshold,upperthreshold);
                obj=disp_hist_xy_sep(obj,30);
            end
        end
        function obj=position_distance_selection(obj,threshold)
            %Can be used before connecting particle positions
            %Finds distances between all particles in each frame
            
            if isempty(obj.pos)
                error('Obtain particle positions first')
            end
            
            input=obj.pos;
            fname=[obj.foldername obj.filename '.tif'];
            info = imfinfo(fname);
            nframes = numel(info);
            pos_new=[];
            
            for n=1:nframes
                
                %Get particle IDs and positions in frame n
                particles=input(find(input(:,3)==n),1:3);
                to_remove=[];
                
                if length(particles)>1
                    %Define empty matrix centers
                    distance=NaN(length(particles(:,1)),length(particles(:,1)));

                    %Fill centers such that each ID has a position in each frame
                    %If the particle was not found that frame, the position is NaN
                    for p=1:length(particles(:,1))
                        for q=1:length(particles(:,1))
                            if p==q
                                distance(p,q)=NaN;
                            else
                                distance(p,q)=sqrt((particles(p,1)-particles(q,1))^2+(particles(p,2)-particles(q,2))^2);
                                if distance(p,q)<threshold
                                    if ~any(to_remove==p)
                                        to_remove=[to_remove ;p];
                                    end
                                    if ~any(to_remove==q)
                                        to_remove=[to_remove ;q];
                                    end
                                end
                            end
                        end
                    end
                    
                    particles_new=[];
                    for p=1:length(particles(:,1))
                        if ~any(to_remove==p)
                            particles_new=[particles_new ; particles(p,:)];
                        end
                    end
                    
                    pos_new=[pos_new ; particles_new];
                    
                end%end if statement n particles larger than 1 in this frame
                
                obj.pos=pos_new;
                
            end%end for loop frames
            
        end
        function obj=remove_particles_flow(obj,lowerthreshold,upperthreshold)
            %function for particles in flow in a microfluidic device
            %removes particles that do not move along with flow (so that
            %are probably stuck)
            
            %Check if object is tracked
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            
            %Find the particles that are stuck
            trace=get_trace(obj);
            nparts=obj.NOP;
            badparticles=[];
            meandisp=[];
            check=10;
            for p=1:nparts
                progress1=p/obj.NOP*100;
                if progress1>check
                    %disp(['Progress in checking for stuck particles: ' num2str(floor(progress1)) '%.'])
                    check=check+10;
                end
                trj=trace(find(trace(:,4)==p),1:3);
                x1=trj(1:length(trj(:,1))-1,1);
                x2=trj(2:length(trj(:,1)),1);
                dx=x2-x1;
                meandisp(p)=nanmean(dx);
                maxdisp(p)=max(dx);
                mindisp(p)=min(dx);
            end
            
            if nargin<2
                meandispall=nanmean(meandisp(:));
                stddispall=nanstd(meandisp(:));
                lowerthreshold=meandispall-nstd*stddispall;
                upperthreshold=meandispall+nstd*stddispall;
            end
            if nargin==2
                disp('Please also provide an upper boundary for the x displacement');
            end
            
            badparticles=find(mindisp<lowerthreshold);
            badparticles=[badparticles find(maxdisp>upperthreshold)];
            
            
            %Remove the particles that are stuck
            trace_old=obj.tr;
            trace_new=[];
            count=0;
            check=10;
            for p=1:obj.NOP
                if ~any(badparticles==p)
                    progress2=p/obj.NOP*100;
                    if progress2>check
                        %disp(['Progress in removing stuck particles: ' num2str(floor(progress2)) '%.'])
                        check=check+10;
                    end
                    
                    count=count+1;
                    trp=trace_old(find(trace_old(:,4)==p),1:3);
                    trace_new=[trace_new ; trp ones(length(trp(:,1)),1)*count];
                end
            end
            %if count~=obj.NOP-length(badparticles)
            %    error('Something went wrong if these values are not the same')
            %end
            obj.NOP=obj.NOP-length(badparticles);
            obj.tr=trace_new;
        end
        function obj=drift_correct(obj,noshow)
            %Takes the trajectories of all particles, finds the average
            %instantanteous displacement, and corrects for drift
            %Noshow determines whether or not before and after correction
            %is shown. (1) is yes, and (0) is no. Default is yes.
            %Drift correction assumes flow on top of brownian motion that
            %is homogenous in time and space
            if nargin<2
                noshow=1;
            end
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            trace=get_trace(obj);
            nframes=obj.NOF;
            nparts=obj.NOP;
            dx_tot=[];
            dy_tot=[];
            %Plot histogram before drift correction
            if noshow==1
                displacement_histogram(obj);
            end
            
            %drift correct
            for p=1:nparts
                trj=trace((p-1)*nframes+1:p*nframes,1:2);
                x1=trj(1:length(trj(:,1))-1,1);
                x2=trj(2:length(trj(:,1)),1);
                dx=x2-x1;
                dx_tot=[dx_tot ; dx];
                y1=trj(1:length(trj(:,1))-1,2);
                y2=trj(2:length(trj(:,1)),2);
                dy=y2-y1;
                dy_tot=[dy_tot ; dy];
            end
            dx_tot=dx_tot(~isnan(dx_tot));
            dy_tot=dy_tot(~isnan(dy_tot));
            dx_correct=mean(dx_tot(:,1));
            dy_correct=mean(dy_tot(:,1));
            for p=1:nparts
                trj=trace((p-1)*nframes+1:p*nframes,1:3);
                trj_size=size(trj);
                trj_new=NaN(trj_size(1),trj_size(2));
                size(trj_new);
                trj_new(:,1)=trj(:,1)-dx_correct*trj(:,3);
                trj_new(:,2)=trj(:,2)-dy_correct*trj(:,3);
                trj_new(:,3)=trj(:,3);
                trace((p-1)*nframes+1:p*nframes,1:3)=trj_new;
            end
            obj.tr=trace;
            
            %Show histogram after drift correct
            if noshow==1
                displacement_histogram(obj);
            end
        end
        function value=find_drift_values(obj,noshow)
            %Finds the average displacement in an object of type particles.
            %If there is constant drift, this displacement divided by the
            %framerate can be considered the drift velocity
            %Outputs value=[driftx drifty]
            
            if nargin<2
                noshow=1;
            end
            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            trace=get_trace(obj);
            nframes=obj.NOF;
            nparts=obj.NOP;
            dx_tot=[];
            dy_tot=[];
            
            %Plot histogram before drift correction
            if noshow==1
                displacement_histogram(obj);
            end
            
            %drift correct
            for p=1:nparts
                trj=trace((p-1)*nframes+1:p*nframes,1:2);
                x1=trj(1:length(trj(:,1))-1,1);
                x2=trj(2:length(trj(:,1)),1);
                dx=x2-x1;
                dx_tot=[dx_tot ; dx];
                y1=trj(1:length(trj(:,1))-1,2);
                y2=trj(2:length(trj(:,1)),2);
                dy=y2-y1;
                dy_tot=[dy_tot ; dy];
            end
            dx_correct=nanmean(dx_tot(:,1));
            dy_correct=nanmean(dy_tot(:,1));
            
            value=[dx_correct dy_correct];
            
        end
        function obj=drift_correct_manual(obj,values)
            %Corrects drift in object based on drift values calculated
            %before
            %values needs to be of form values=[dx dy] where dx is the
            %drift displacement in x over one frame and idem dito for dy.
            %WARNING: overrides obj.tr. Do not run this function if you
            %are not sure about the driftvalues or if you have no backup of
            %the trace.
            
            if nargin<2
                disp('This function only works if you provide driftvalues.')
                error('If you want to do automatic drift correct, try "drift_correct" instead.')
            end
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            nparts=obj.NOP;
            nframes=obj.NOF;
            trace=obj.tr;
            
            dx_correct=values(1);
            dy_correct=values(2);
            
            for p=1:nparts
                trj=trace((p-1)*nframes+1:p*nframes,1:3);
                trj_size=size(trj);
                trj_new=NaN(trj_size(1),trj_size(2));
                size(trj_new);
                trj_new(:,1)=trj(:,1)-dx_correct*trj(:,3);
                trj_new(:,2)=trj(:,2)-dy_correct*trj(:,3);
                trj_new(:,3)=trj(:,3);
                trace((p-1)*nframes+1:p*nframes,1:3)=trj_new;
            end
            obj.tr=trace;
            drift_corr_hist_add=[dx_correct,dy_correct];
            obj.drift_corr_hist=[obj.drift_corr_hist drift_corr_hist_add'];
            
        end
        function obj=make_dist_matrix(obj)
            %Takes the positions in each frame and finds a distance matrix dist
            %dist is n x m x m where m is number of particles and n is
            %number of frames
            %Uses trace as an input, so define obj.tr first
            %Note that for long movies with many particles this script is
            %slow
            
            %if trace does not exist, trace particles first
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            
            input=obj.tr;
            nparts=obj.NOP;
            nframes=obj.NOF;
            distance_all=NaN(obj.NOP,obj.NOP,obj.NOF);
            %distance_all=cell(1,1);
            count=0;
            if max(input(:,4))<2
                obj.dist=distance_all;
                disp('You only have one particle. There are no distances to be calculated');
            else
                
                
                for n=1:nframes
                    
                    %Get particle IDs and positions in frame n
                    particles=[input(find(input(:,3)==n),1:2) input(find(input(:,3)==n),4)];
                    
                    %Let user know how far the code is
                    if ~mod(n/100,1)
                        disp(['Working on frame ' num2str(n) ' out of ' num2str(obj.NOF)])
                    end
                    
                    %Define empty matrix centers
                    centers=NaN(nparts,3);
                    
                    %Fill centers such that each ID has a position in each frame
                    %If the particle was not found that frame, the position is NaN
                    for p=1:nparts
                        if ismember(p,particles(:,3))
                            centers(p,1:3)=particles(find(particles(:,3)==p),1:3);
                        else
                            centers(p,3)=p;
                        end
                    end
                    
                    nparts=max(input(:,4));
                    distance=NaN(nparts,nparts);
                    for i=1:length(centers(:,1))
                        for k=1:length(centers(:,1))
                            %Get distance matrix
                            distance(centers(i,3),centers(k,3))=sqrt((centers(i,1)-centers(k,1))^2+(centers(i,2)-centers(k,2))^2);
                        end
                    end
                    count=count+1;
                    distance_all(:,:,n)=distance;
                    %distance_all{count}=distance;
                end
                obj.dist=distance_all;
            end
        end
        function obj=make_av_dist_matrix(obj)
            %Calculates the average distance matrix
            
            %if trace does not exist, trace particles first
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            
            input=obj.tr;
            nparts=obj.NOP;
            nframes=obj.NOF;
            distance_av_old=zeros(obj.NOP,obj.NOP);
            count=0;
            if max(input(:,4))<2
                obj.dist=distance_all;
                disp('You only have one particle. There are no distances to be calculated');
            else
                
                
                for n=1:nframes
                    
                    %Get particle IDs and positions in frame n
                    particles=[input(find(input(:,3)==n),1:2) input(find(input(:,3)==n),4)];
                    
                    %Let user know how far the code is
                    if ~mod(n/100,1)
                        disp(['Working on frame ' num2str(n) ' out of ' num2str(obj.NOF)])
                    end
                    
                    %Define empty matrix centers
                    centers=NaN(nparts,3);
                    
                    %Fill centers such that each ID has a position in each frame
                    %If the particle was not found that frame, the position is NaN
                    for p=1:nparts
                        if ismember(p,particles(:,3))
                            centers(p,1:3)=particles(find(particles(:,3)==p),1:3);
                        else
                            centers(p,3)=p;
                        end
                    end
                    
                    nparts=max(input(:,4));
                    distance=NaN(nparts,nparts);
                    for i=1:length(centers(:,1))
                        for k=1:length(centers(:,1))
                            %Get disttance matrix
                            distance(centers(i,3),centers(k,3))=sqrt((centers(i,1)-centers(k,1))^2+(centers(i,2)-centers(k,2))^2);
                        end
                    end
                    %if ~isnan(distance)
                    distance_av=(count*distance_av_old+distance)/(count+1);
                    %distance_av(isnan(distance_av))=distance_av_old(isnan(distance_av));
                    distance_av_old(find(~isnan(distance_av)))=distance_av(~isnan(distance_av));
                    %distance_av_old=distance_av;
                    count=count+1;
                    %end
                end
                try obj.dist_av=distance_av_old;
                catch
                    obj.dist_av=NaN;
                    disp('For some reason all distances were NaN. Check make_dist_av code')
                end
            end
        end
        function plot_trace_on_movie(obj,start_frame, end_frame)
            %A is the total stack of images corresponding to the track
            %Displays a movie of the droplets in each frame with the
            %corresponding position on top.
            %Settings for video creation
            %Saves the video to the indicated foldername
            
            if isnan(obj.im)
                error('There is no imagestack loaded in the object, so no trajectory can be plotted')
            end
            %if trace does not exist, trace particles first
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
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
                if tmax>nframes
                    tmax=nframes;
                end
                
                for k=t+1:tmax-1
                    %read image
                    fname=[obj.foldername obj.filename '.tif'];
                    Ik=imread(fname,k);
                    
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
                    %plot flow profile
                    for i=1:nparts
                        quiver(x2(k,i),y2(k,i),factor*u2(k,i),factor*v2(k,i),linesize,'Color','r','LineWidth',linewidth)
                    end
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
        function bad_frames=find_abs_displacement(obj,noshow)
            %calculates the absolute displacement per particle in each frame and
            %returns:
            % - a scattereplot with a different color for each particle, containing the absolute
            %   displacement versus the framenumber
            % - a list of bad frames
            %   bad frames are defined as frames where the displacement is
            %   larger than the average + 3 standard deviations.
            % the input parameter noshow indicates whether or not to show
            % the graphs. Input 1 if yes, and 0 if no.
            if nargin<2
                noshow=0;
            end
            
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            trace=get_trace(obj);
            nparts=obj.NOP;
            bad_frames=cell(nparts,1);
            cmap=jet(nparts);
            if noshow==1
                figure
                hold on
                box on
                lgd=legend('show');
            end
            for p=1:nparts
                trj=trace(find(trace(:,4)==p),1:3);
                x1=trj(1:length(trj(:,1))-1,1);
                x2=trj(2:length(trj(:,1)),1);
                y1=trj(1:length(trj(:,1))-1,2);
                y2=trj(2:length(trj(:,1)),2);
                dx=(x2-x1).^2;
                dy=(y2-y1).^2;
                abs_disp=sqrt(dx+dy);
                if noshow==1
                    plot(trj(1:length(trj(:,1))-1,3), abs_disp, 'o','Color',cmap(p,:),'DisplayName',['Particle ' num2str(p)])
                end
                averagedisp=mean(abs_disp);
                stdevdisp=std(abs_disp);
                bad_frames{p}=find(abs_disp>averagedisp+3*stdevdisp);
            end
            if noshow==1
                lgd.Color=[0.5 0.5 0.5];
                lgd.TextColor=[1 1 1];
                ylabel('Absolute displacement (pixels)','FontSize',18)
                xlabel('Framenumber','FontSize',18)
                hold off
            end
        end
        function obj=complete_trace(obj)
            %Finds holes in trace and fills them with NaN
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            trace=obj.tr;
            nparts=obj.NOP;
            nframes=obj.NOF;
            trace_new=[];
            for p=1:nparts
                trj=trace(find(trace(:,4)==p),1:4);
                count=0;
                if length(trj(:,1))<nframes
                    trj_new=NaN(nframes,4);
                    for n=1:nframes
                        if any(trj(:,3)==n)
                            trj_new(n,1:4)=trj(n-count,1:4);
                        else
                            trj_new(n,1:4)=[NaN NaN n p];
                            count=count+1;
                        end
                    end
                    trace_new=[trace_new ; trj_new];
                else
                    trace_new=[trace_new ; trj];
                end
                
            end
            obj.tr=trace_new;
        end
        function obj=automatically_improve_trace(obj)
            %takes the trace contained in object, finds frames where the
            %displacement is anomalously large, removes those frames, and
            %tracks again. This removes errors due to incorrectly finding
            %particles.
            if isempty(obj.tr)
                obj=track(obj);
                disp('Executed tracking first because tr was still empty')
            end
            
            obj=complete_trace(obj);
            
            bad_frames=find_abs_displacement(obj,0);
            nparts=obj.NOP;
            nframes=obj.NOF;
            trace=obj.tr;
            tr_temp=NaN(nparts*nframes,4);
            for p=1:nparts
                trj=trace(find(trace(:,4)==p),1:4);
                badframesp=bad_frames{p};
                trj(badframesp,1:2)=NaN;
                tr_temp(p:nparts:(nparts)*nframes-nparts+p,:)=trj(:,:);
            end
            
            maxdisp=get_numerical_input('What is the max displacement? ');
            param.mem=get_numerical_input('How long can a particle go missing? Choose value larger than 1: ');
            param.dim=2; % dimensionality of data
            param.good=get_numerical_input('What is the minimal number of frames for a track to be accepted? ');
            param.quiet=1; % 0 = text, 1 = no text
            obj.tr=track(tr_temp(:,1:3),maxdisp,param);
            obj=complete_trace(obj);
            
        end
        function plot_tracked_image(obj,frame,vectorsize,linewidth)
            %Plots framenumber 'frame' of the associated tiffstack and
            %overlays the found particles in this image
            %Default frame is 1;
            %Default vectorsize is 3;
            %Default linewidth is 2;
            
            if nargin==1
                frame=1;
                vectorsize=1;
                linewidth=2;
            elseif nargin==2
                linewidth=2;
                vectorsize=1;
            elseif nargin==3
                linewidth=2;
            end
            
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            trace=obj.tr;
            
            fname=[obj.foldername obj.filename '.tif'];
            I=imread(fname,frame);
            x2=NaN(1,obj.NOP);
            y2=NaN(1,obj.NOP);
            u2=NaN(1,obj.NOP);
            v2=NaN(1,obj.NOP);
            
            for i=1:obj.NOP
                trj=trace(find(trace(:,4)==i),1:3);
                x2(1,i)=trj(frame,1);
                y2(1,i)=trj(frame,2);
                %calculate x component motion particle
                u2(1,i)=trj(frame+1,1)-trj(frame,1);
                %calculate y component motion particle
                v2(1,i)=trj(frame+1,2)-trj(frame,2);
            end
            
            figure
            hold on
            box on
            %show original image in grayscale
            colormap gray
            imagesc(I);
            sizeI=size(I);
            %show vectorplot of droplet motion in red
            for i=1:obj.NOP
                quiver(x2(1,i),y2(1,i),u2(1,i),v2(1,i),vectorsize,'Color','r','LineWidth',linewidth,'MaxHeadSize',1)
            end
            %Settings aspect ratios
            set(gca,'Ylim',[0 sizeI(1)],'Xlim', [0 sizeI(2)])
            %set(gca,'Ytick',[0:200:sizeI(2)],'Xtick',[0:300:sizeI(1)])
            pbaspect([sizeI(1)/sizeI(1) sizeI(1)/sizeI(2) 1]);
            %pbaspect([1 1 1])
            %axis([780 1280 300 800])
            set(gca,'YDir','normal');
            
        end
        function obj=remove_point_trace(obj,frame)
            % manually change the information in trace to NaN for one frame
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            trace=obj.tr;
            trace(find(trace(:,3)==frame),1:2)=NaN;
            obj.tr=trace;
        end
        function obj=calc_VCF(obj,n)
            %velocity correlation function
            %n is the fraction of the trace used as the largest time
            %interval. Default n is 3
            %created on 17-06-15 by Pepijn Moerman
            %last modified on 04-07 to fit in tracking object
            %purpose: calculate the correlation in angle
            %between steps in a trajectory of active particles
            %as function of the time lag between steps
            
            if nargin==1
                n=3;
            end
            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            nof=obj.NOF;
            trace=obj.tr;
            l=floor(nof/n);
            phi = NaN(l+1,obj.NOP+1);
            phi(:,1)=[0:l]'/obj.framerate;
            
            %loop over all particles
            for p=1:obj.NOP
                %collect particle trajectory and define matrices
                trj=trace(find(trace(:,4)==p),1:3);
                
                %loop over all time intervals
                for d=0:l
                    theta=NaN(nof,1);
                    
                    for i=1:d:nof-2-(d+1)
                        x1= trj(i,1);
                        y1= trj(i,2);
                        x2 = trj(i+1,1);
                        y2 = trj(i+1,2);
                        x3 = trj(i+d,1);
                        y3 = trj(i+d,2);
                        x4 = trj(i+d+1,1);
                        y4 = trj(i+d+1,2);
                        
                        xn1=x2-x1;
                        yn1=y2-y1;
                        xn2=x4-x3;
                        yn2=y4-y3;
                        
                        a=(yn2-yn1)^2+(xn2-xn1)^2;
                        b=xn1^2+yn1^2;
                        c=xn2^2+yn2^2;
                        
                        theta(i,1) = (b+c-a)/(2*sqrt(b*c));
                    end % end loop over time for average
                    
                    
                    phi(d+1,p+1)=nanmean(theta(:,1));
                end %end loop over time interval
                
            end %end loop over particle number
            
            obj.VCF=phi;
            
            figure
            hold on
            box on
            cmap=jet(obj.NOP);
            for p=1:obj.NOP
                %plot(phi(:,1),phi(:,p+1),'.-','MarkerSize',10,'Color',cmap(p,:))
                plot(phi(:,1),phi(:,p+1),'k.-','MarkerSize',10)
            end
            axis([0 l/obj.framerate -1.1, 1.2])
            xlabel('\delta t (s)','FontSize',14)
            ylabel('C (\deltat)','FontSize',14)
            hold off
        end
        function obj=calc_sint(obj,n)
            %calculates a variation on the velocity correlation function
            %defined as 3sin(theta)^2-2
            %this function is 0 if the angle is randomly distributed and 1
            %if the angle is perpendicular and -1 for a parallel angle.
            %As opposed to VCF, which is 1 for parallel angle and 0 for
            %either a random or perpendicular angle.
            %n is the fraction of the trace used as the largest time
            %interval. Default n is 3
            %created on 17-06-15 by Pepijn Moerman
            %last modified on 04-07 to fit in tracking object
            %purpose: calculate the correlation in angle
            %between steps in a trajectory of active particles
            %as function of the time lag between steps
            
            if nargin==1
                n=3;
            end
            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            nof=obj.NOF;
            trace=obj.tr;
            l=round(nof/n);
            phi = NaN(l+1,obj.NOP+1);
            phi(:,1)=[0:l]'/obj.framerate;
            
            %loop over all particles
            for p=1:obj.NOP
                %collect particle trajectory and define matrices
                trj=trace(find(trace(:,4)==p),1:3);
                
                %loop over all time intervals
                for d=0:l
                    theta=NaN(nof,1);
                    
                    for i=1:nof-2-(d+1)
                        x1= trj(i,1);
                        y1= trj(i,2);
                        x2 = trj(i+1,1);
                        y2 = trj(i+1,2);
                        x3 = trj(i+d,1);
                        y3 = trj(i+d,2);
                        x4 = trj(i+d+1,1);
                        y4 = trj(i+d+1,2);
                        
                        xn1=x2-x1;
                        yn1=y2-y1;
                        xn2=x4-x3;
                        yn2=y4-y3;
                        
                        a=(yn2-yn1)^2+(xn2-xn1)^2;
                        b=xn1^2+yn1^2;
                        c=xn2^2+yn2^2;
                        
                        theta(i,2) = (b+c-a)/(2*sqrt(b*c));
                    end % end loop over time for average
                    
                    theta_nonan=theta(:,2);
                    theta_nonan(isnan(theta_nonan(:,1)),:)=[];
                    
                    phi(d+1,p+1)=mean(theta_nonan(:,1));
                end %end loop over time interval
                
            end %end loop over particle number
            
            obj.VCF=phi;
            
            figure
            hold on
            box on
            cmap=jet(obj.NOP);
            for p=1:obj.NOP
                plot(phi(:,1),phi(:,p+1),'.-','MarkerSize',10,'Color',cmap(p,:))
            end
            axis([0 l/obj.framerate -1.1, 1.2])
            xlabel('time (s)','FontSize',16)
            ylabel('C (\deltat)','FontSize',18)
            hold off
        end
        function obj=average_VCF(obj)
            if isempty(obj.VCF)
                obj=calc_VCF(obj);
            end
            
            if obj.NOP>1
                average_VCF=nanmean(obj.VCF(:,2:obj.NOP),2);
            else
                average_VCF=obj.VCF(:,2);
            end
            obj.VCFav=average_VCF;
            
            figure
            hold on
            box on
            plot(obj.VCF(:,1)/obj.framerate,average_VCF(:,1),'.-','MarkerSize',20,'Color','Black')
            axis([0 obj.VCF(end,1) -1.1, 1.2])
            xlabel('time (s)','FontSize',16)
            ylabel('C (\deltat)','FontSize',18)
            hold off
            
        end
        function obj=merge_tracks(obj,ID1,ID2)
            %manually merges tracks of particles with ID labels ID1 and ID2
            %Only merges trajectories with no overlapping framenumbers
            %plots trajectories before and after merger
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            trace=obj.tr;
            plot_tracks(obj);
            further=1;
            trj1=trace(find(trace(:,4)==ID1),1:3);
            trj1=trj1(find(~isnan(trj1(:,1))),1:3);
            trj2=trace(find(trace(:,4)==ID2),1:3);
            trj2=trj2(find(~isnan(trj2(:,1))),1:3);
            if max(trj1(:,3))<min(trj2(:,3))
                trj=[trj1(:,1:3) ; trj2(:,1:3)];
            elseif min(trj1(:,3))>max(trj2(:,3))
                trj=[trj2(:,1:3) ; trj1(:,1:3)];
                sizies=size(trj);
            else
                obj.tr=trace;
                further=0;
            end
            if further~=0
                done=1;
                count=0;
                trace_new=[];
                for p=1:obj.NOP
                    count=count+1;
                    if any([ID1 ID2]==p)
                        if done==1
                            trace_new=[trace_new ; trj(:,1:3) ones(length(trj(:,1)),1)*count];
                            done=0;
                        end
                    else
                        trace_new=[trace_new ; trace(find(trace(:,4)==p),1:4)];
                    end
                end
                obj.tr=trace_new;
                obj.NOP=max(trace_new(:,4));
                obj=complete_trace(obj);
                plot_tracks(obj);
            end
        end
        function obj=calc_MSD(obj,fraction,selection)
            %calculates mean squared displacement per particle and the
            %average mean squared displacement.
            %fraction indicated what fraction of the trajectory to use to
            %calculated MSD (last few steps have poor statistics). Default
            %is 1/3; i.e. fraction is 10;
            %selection indicates for which particles to do calculation
            %format must be [ID1 ID2 ID3...]
            %if selection is not provided or zero, all particles are used in calculation.
            
            
            %!!!IMPORTANT NOTICE: at 9-10-18 this code was changed to save
            %data in frames and pixels instead of um and s for reasons of
            %unity in the object. All instances of this object created
            %before that date need to be recalculated to get correct MSD
            %values!!!
            
            if nargin<2
                fraction=10;
            end
            if nargin<3
                selection=1:obj.NOP;
            end
            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            trace=obj.tr;
            nparts=obj.NOP;
            nframes=obj.NOF;
            framerate_temp=obj.framerate;
            scale_temp=obj.scale;
            obj.MSD=cell(nparts+1,2);
            
            if length(trace(:,1))~=max(trace(:,4))*max(trace(:,3))
                obj=complete_trace(obj);
            end
            
            MSD_tot=NaN(floor(nframes/fraction)+1,nparts);
            for i=1:nparts
                % Select the part of the file about particle i
                trj = trace((i-1)*nframes+1:i*nframes,1:3);
                if ~isempty(trj)
                    
                    % Find the length of the trajectory of particle i
                    tmax=length(trj(:,1));
                    
                    % Define empty matrix for MSD
                    MSD_temp=NaN(floor(tmax/fraction)+1,nparts);
                    MSD_temp(1,1:2)=[0 0];
                    
                    % For a time interval of 1 to 1/3 of the length of the trajectory
                    for t=1:floor(tmax/fraction)
                        %Create two matrixes with a difference t between them
                        a_frontstep=trj(t+1:t:length(trj),1:2);
                        a_backstep=trj(1:t:length(trj)-t,1:2);
                        %Substract to obtain MSD
                        MSD_all=(a_frontstep(:,1)-a_backstep(:,1)).^2+(a_frontstep(:,2)-a_backstep(:,2)).^2;
                        %Average to obtain average MSD for t and save
                        MSD_temp(t+1,1)=t;
                        MSD_temp(t+1,2)=nanmean(MSD_all(:,1));
                    end
                    
                    %End if statements if matrix is empty after removing stuck particles
                end
                %End loop over particle number
                obj.MSD{1+i,1}=i;
                obj.MSD{1+i,2}=MSD_temp(:,1:2);
                MSD_tot(:,i)=MSD_temp(:,2);
            end
            obj.MSD{1,1}='Average';
            obj.MSD{1,2}=[MSD_temp(:,1) mean(MSD_tot(:,:),2)];
            
        end
        function plot_MSD(obj,fraction)
            %plots the Mean Squared Displacement for each particle in time
            %plots the average Mean Squared Displacement
            %fraction indicated what fraction of the trajectory to use to
            %calculated MSD (last few steps have poor statistics). Default
            %is 1/3; i.e. fraction is 3;
            
            %!!!IMPORTANT NOTICE: at 9-10-18 this code was changed to save
            %data in frames and pixels instead of um and s for reasons of
            %unity in the object. All instances of this object created
            %before that date need to be recalculated to get correct MSD
            %values!!!
            
            if nargin<2
                fraction=3;
            end
            
            if isempty(obj.MSD) || fraction~=3
                obj=calc_MSD(obj,fraction);
            end
            if obj.NOP>1
                %Plot the individual MSDs in one figure
                figure
                cmap=jet(obj.NOP);
                hold on
                box on
                
                for i=1:obj.NOP
                    MSDi=obj.MSD{i+1,2};
                    particleID=obj.MSD{i+1,1};
                    plot(MSDi(:,1)/obj.framerate,MSDi(:,2)*obj.scale^2,'o','Color',cmap(i,:),'DisplayName',['Particle ' num2str(particleID)])
                end
                xlabel('\Deltat (s)','FontSize',18)
                ylabel('MSD (\mum^2/s)','FontSize',18)
                lgd=legend('show');
                lgd.Color=[0.5 0.5 0.5];
                lgd.TextColor=[1 1 1];
                hold off
            end
            
            %Plot the average MSD in one figure
            figure
            MSD_average=obj.MSD{1,2};
            hold on
            box on
            %lgd=legend('show');
            %plot(MSD_average(:,1),MSD_average(:,2),'ko','DisplayName','Average')
            plot(MSD_average(:,1)/obj.framerate,MSD_average(:,2)*obj.scale^2,'ko')
            xlabel('\Deltat (s)','FontSize',18)
            ylabel('MSD (\mum^2/s)','FontSize',18)
            %lgd.Color=[0.5 0.5 0.5];
            %lgd.TextColor=[1 1 1];
            hold off
            
            %Plot the aveage MSD in log log representation
            figure
            MSD_average=obj.MSD{1,2};
            box on
            loglog(MSD_average(:,1)/obj.framerate,MSD_average(:,2)*obj.scale^2,'ko')
            xlabel('\Deltat (s)','FontSize',18)
            ylabel('MSD (\mum^2/s)','FontSize',18)
        end
        function badparticles=find_stagnant_particles(obj,n)
            %Takes a set of tracked particles and finds the sum of the
            %absolute displacement. If the sum is lower than 3*stdev, it is
            %considered stuck.
            %n indicates how many times the stdev the displacement needs to
            %be below the average displacement to be stuck.
            %default is 1
            if nargin<2
                n=3;
            end
            
            
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
            trace=get_trace(obj);
            nparts=obj.NOP;
            badparticles=[];
            sumdisp=[];
            
            for p=1:nparts
                trj=trace(find(trace(:,4)==p),1:3);
                x1=trj(1:length(trj(:,1))-1,1);
                x2=trj(2:length(trj(:,1)),1);
                y1=trj(1:length(trj(:,1))-1,2);
                y2=trj(2:length(trj(:,1)),2);
                dx=(x2-x1).^2;
                dy=(y2-y1).^2;
                abs_disp=sqrt(dx+dy);
                
                sumdisp(p)=nansum(abs_disp);
            end
            avdisp=mean(sumdisp);
            stdisp=std(sumdisp);
            badparticles=find(sumdisp<avdisp-n*stdisp);
            
        end
        function obj=remove_stagnant_particles(obj,n)
            %Takes a tracked object and removes stagnant particles from it
            
            if nargin<2
                n=3;
            end
            badparticles=find_stagnant_particles(obj,n);
            
            trace_old=obj.tr;
            trace_new=[];
            count=0;
            for p=1:obj.NOP
                if ~any(badparticles==p)
                    count=count+1;
                    trp=trace_old(find(trace_old(:,4)==p),1:3);
                    trace_new=[trace_new ; trp ones(length(trp(:,1)),1)*count];
                end
            end
            if count~=obj.NOP-length(badparticles)
                error('Something went wrong if these values are not the same')
            end
            obj.NOP=obj.NOP-length(badparticles);
            obj.tr=trace_new;
            
        end
        function obj=find_reduced_trajectory(obj,dt,noshow)
            %This function gives the tracks around the frame of reference
            %of particle each central particle
            %The input variable noshow indicates whether or not to plot the
            %reduced tracks. Give (1) for yes and (0) for no. Default is
            %no.
            %dt is the stepsize over which you want to measure the reduced
            %trajectory. Default is 1
            if nargin<3
                noshow=0;
            end
            if nargin<2
                dt=1;
            end
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            trace=obj.tr;
            nof=obj.NOF;
            if nof<3
                disp('There are fewer than 3 frames. That is not enough to determine the rotation angle.')
                obj.tr_red=NaN;
            else
                obj.tr_red={};
                %For each particle
                for p=1:obj.NOP
                    %Define three trajectories with a difference of 1
                    %timestep
                    trace_p=trace(find(trace(:,4)==p),1:2);
                    trj_0=trace_p(1+dt:dt:nof-dt,1:2);
                    trj_1=trace_p(1:dt:nof-2*dt,1:2);
                    trj_2=trace_p(1+2*dt:dt:nof,1:2);
                    
                    %Define the reduced trajectories
                    trj_0_red=zeros(length(trj_0(:,1)),2);
                    trj_1_red=trj_1(:,1:2)-trj_0(:,1:2);
                    trj_2_red=trj_2(:,1:2)-trj_0(:,1:2);
                    
                    %Rotate the reduced trajectories to align with the
                    %x-axis
                    trj_1_rot=[];
                    trj_2_rot=[];
                    for j=1:length(trj_0_red)
                        angle=pi-atan2(trj_1_red(j,2),trj_1_red(j,1));
                        R=[cos(angle)  -sin(angle) ; sin(angle)  cos(angle)];
                        trj_1_rot(j,1:2)=R*trj_1_red(j,1:2)';
                        trj_2_rot(j,1:2)=R*trj_2_red(j,1:2)';
                    end
                    if noshow==1
                        figure
                        hold on
                        box on
                        axis equal
                        xlabel('X (\mum)')
                        ylabel('Y (\mum)')
                        plot(trj_0_red(:,1)*obj.scale,trj_0_red(:,2)*obj.scale,'ko')
                        plot(trj_1_rot(:,1)*obj.scale,trj_1_rot(:,2)*obj.scale,'r-')
                        plot(trj_2_rot(:,1)*obj.scale,trj_2_rot(:,2)*obj.scale,'b-')
                        hold off
                    end
                    tr_red_temp=[trj_0_red(:,1:2) ones(length(trj_0_red(:,1)),1)*0 ; trj_1_rot(:,1:2) ones(length(trj_1_rot(:,1)),1)*1 ; trj_2_rot(:,1:2) ones(length(trj_2_rot(:,1)),1)*2];
                    obj.tr_red{p}=tr_red_temp;
                    
                end
            end
        end
        function obj=plot_reduced_tracks(obj,dt,cutoffvalue)
            %plots the reducted trajectory of each particle in the chain
            %with 2 neighbours
            %Cutoffvalue gives distance between points used to give
            %determine color scheme
            %Default cutoffvalue=2;
            %dt is the stepsize for which reduced tracks are calculated
            %default is 1;
            if nargin==1
                dt=1;
                cutoffvalue=2;
            elseif nargin==2
                cutoffvalue=2;
            end
            
            
            obj=find_reduced_trajectory(obj,dt,0);
            
            %For each particle
            for p=1:obj.NOP
                %Obtain the reduced trajectories
                tr_red_p=obj.tr_red{p};
                trj_0_red=tr_red_p(find(tr_red_p(:,3)==0),1:2);
                trj_1_red=tr_red_p(find(tr_red_p(:,3)==1),1:2);
                trj_2_red=tr_red_p(find(tr_red_p(:,3)==2),1:3);
                
                %Get color scheme based on number of neighbours for
                %ID2;
                points=chains(get_name(obj),get_parameters(obj),'n');
                points=add_pos(points,trj_2_red(:,1:3));
                points_trace=[trj_2_red(:,1:2) ones(length(trj_2_red(:,1)),1) ones(length(trj_2_red(:,1)),1)];
                for i=1:length(points_trace(:,1))
                    points_trace(i,4)=i;
                end
                points=add_trace(points,points_trace);
                points=find_neighbours(points,cutoffvalue);
                num_of_neighs=points.num_neighs;
                max_num=max(num_of_neighs);
                if max_num<1
                    max_num=1;
                end
                cmap=jet(max_num);
                
                %Plot the reduced trajectories
                figure
                hold on
                box on
                axis equal
                xlabel('X (\mum)')
                ylabel('Y (\mum)')
                title(['Particle ' num2str(p)])
                plot(trj_0_red(:,1)*obj.scale,trj_0_red(:,2)*obj.scale,'ko')
                %plot(trj_1_red(:,1)*obj.scale,trj_1_red(:,2)*obj.scale,'ro')
                %plot the last reduced trajectory in colors
                %depending on the number of neighbours
                for i=1:length(trj_2_red(:,1))
                    num_neigh_i=round(num_of_neighs(i));
                    if num_neigh_i<1
                        num_neigh_i=1;
                    elseif num_neigh_i>max_num
                        num_neigh_i=max_num;
                    end
                    plot(trj_2_red(i,1)*obj.scale,trj_2_red(i,2)*obj.scale,'o','Color',cmap(num_neigh_i,:))
                    %plot(trj_2_red(i,1)*obj.scale,trj_2_red(i,2)*obj.scale,'o')
                end
                hold off
            end
            
        end
        function output=plot_average_reduced_track(obj,dt,cutoffvalue,axislimits)
            %plots the reducted trajectory of each particle in the chain
            %with 2 neighbours
            %Cutoffvalue gives distance between points used to give
            %determine color scheme
            %Default cutoffvalue=2;
            if nargin==1
                cutoffvalue=2;
                dt=1;
            elseif nargin==2
                cutoffvalue=2;
            end
            
            obj=find_reduced_trajectory(obj,dt,0);
            
            trj_0_red=[];
            trj_1_red=[];
            trj_2_red=[];
            for p=1:obj.NOP
                tr_red_p=obj.tr_red{p};
                %find trajectories of particle p
                trj_0_red_p=tr_red_p(find(tr_red_p(:,3)==0),1:2);
                trj_1_red_p=tr_red_p(find(tr_red_p(:,3)==1),1:2);
                trj_2_red_p=tr_red_p(find(tr_red_p(:,3)==2),1:3);
                
                %append into trajectories of all particle
                trj_0_red=[trj_0_red ; trj_0_red_p];
                trj_1_red=[trj_1_red ; trj_1_red_p];
                trj_2_red=[trj_2_red ; trj_2_red_p];
            end
            %Get color scheme based on number of neighbours for
            %ID2;
            points=chains(get_name(obj),get_parameters(obj),'n');
            points=add_pos(points,trj_2_red(:,1:3));
            points_trace=[trj_2_red(:,1:2) ones(length(trj_2_red(:,1)),1) ones(length(trj_2_red(:,1)),1)];
            %points_trace=points_trace(~isnan(points_trace(:,1)),:);
            for i=1:length(points_trace(:,1))
                points_trace(i,4)=i;
            end
            points=add_trace(points,points_trace);
            points=make_av_dist_matrix(points);
            points=find_neighbours(points,cutoffvalue);
            num_of_neighs=points.num_neighs;
            max_num=max(num_of_neighs);
            cmap=jet(max_num);
            output=points;
            if nargin<4
                if min(trj_2_red(:,1))<0
                    xmin=1.3*min(trj_2_red(:,1));
                else
                    xmin=0.7*min(trj_2_red(:,1));
                end
                if min(trj_2_red(:,2))<0
                    ymin=1.3*min(trj_2_red(:,2));
                else
                    ymin=0.7*min(trj_2_red(:,2));
                end
                xmax=1.3*max(trj_2_red(:,1));
                ymax=1.3*max(trj_2_red(:,2));
                limits=max([abs(xmin) abs(xmax) abs(ymin) abs(ymax)]);
                %axislimits=[xmin 1.3*max(trj_2_red(:,1)) ymin 1.3*max(trj_2_red(:,2))]*obj.scale;
                axislimits=[-limits limits -limits limits]*obj.scale;
            end
            
            %Plot the reduced trajectories
            figure
            hold on
            box on
            axis equal
            axis(axislimits)
            xlabel('X (\mum)')
            ylabel('Y (\mum)')
            title(['Particle ' num2str(p)])
            %plot(trj_0_red(:,1)*obj.scale,trj_0_red(:,2)*obj.scale,'ko')
            %plot(trj_1_red(:,1)*obj.scale,trj_1_red(:,2)*obj.scale,'ro')
            %plot the last reduced trajectory in colors
            %depending on the number of neighbours
            for i=1:length(trj_2_red(:,1))
                num_neigh_i=round(num_of_neighs(i));
                if num_neigh_i<1
                    num_neigh_i=1;
                end
                plot(trj_2_red(i,1)*obj.scale,trj_2_red(i,2)*obj.scale,'o','Color',cmap(num_neigh_i,:))
            end
            line([0 0],[-1000 1000],'Color','Black');
            line([-1000 1000],[0 0],'Color','Black');
            hold off
            
        end
        function obj=plot_average_reduced_track_subplot(obj,dt,cutoffvalue,axislimits)
            %plots the reducted trajectory of each particle in the chain
            %with 2 neighbours
            %Cutoffvalue gives distance between points used to give
            %determine color scheme
            %Default cutoffvalue=2;
            %dt must be a list of times ordered as [t1 t2 t3 t4]
            if nargin==1
                cutoffvalue=2;
                dt=1;
            elseif nargin==2
                cutoffvalue=2;
            end
            
            figure
            box on
            
            for q=1:length(dt)
                disp(['Starting with time ' num2str(q) ' of ' num2str(length(dt)) '.'])
                dti=dt(q);
                output=find_reduced_trajectory_temp(obj,dti,0);
                
                trj_0_red=[];
                trj_1_red=[];
                trj_2_red=[];
                for p=1:obj.NOP
                    tr_red_p=output{p};
                    %find trajectories of particle p
                    trj_0_red_p=tr_red_p(find(tr_red_p(:,3)==0),1:2);
                    trj_1_red_p=tr_red_p(find(tr_red_p(:,3)==1),1:2);
                    trj_2_red_p=tr_red_p(find(tr_red_p(:,3)==2),1:3);
                    
                    %append into trajectories of all particle
                    trj_0_red=[trj_0_red ; trj_0_red_p];
                    trj_1_red=[trj_1_red ; trj_1_red_p];
                    trj_2_red=[trj_2_red ; trj_2_red_p];
                end
                %Get color scheme based on number of neighbours for
                %ID2;
                points=chains(get_name(obj),get_parameters(obj),'n');
                points=add_pos(points,trj_2_red(:,1:3));
                points_trace=[trj_2_red(:,1:2) ones(length(trj_2_red(:,1)),1) ones(length(trj_2_red(:,1)),1)];
                for j=1:length(points_trace(:,1))
                    points_trace(j,4)=j;
                end
                points=add_trace(points,points_trace);
                points=find_neighbours(points,cutoffvalue);
                num_of_neighs=points.num_neighs;
                max_num=max(num_of_neighs);
                cmap=jet(max_num);
                
                if nargin<4
                    if min(trj_2_red(:,1))<0
                        xmin=1.3*min(trj_2_red(:,1));
                    else
                        xmin=0.7*min(trj_2_red(:,1));
                    end
                    if min(trj_2_red(:,2))<0
                        ymin=1.3*min(trj_2_red(:,2));
                    else
                        ymin=0.7*min(trj_2_red(:,2));
                    end
                    xmax=1.3*max(trj_2_red(:,1));
                    ymax=1.3*max(trj_2_red(:,2));
                    limits=max([abs(xmin) abs(xmax) abs(ymin) abs(ymax)]);
                    %axislimits=[xmin 1.3*max(trj_2_red(:,1)) ymin 1.3*max(trj_2_red(:,2))]*obj.scale;
                    axislimits=[-limits limits -limits limits]*obj.scale;
                end
                
                %Plot the reduced trajectories
                subplot(1,length(dt),q)
                %plot(trj_0_red(:,1)*obj.scale,trj_0_red(:,2)*obj.scale,'ko')
                %plot(trj_1_red(:,1)*obj.scale,trj_1_red(:,2)*obj.scale,'ro')
                %plot the last reduced trajectory in colors
                %depending on the number of neighbours
                for j=1:length(trj_2_red(:,1))
                    subplot(1,length(dt),q)
                    hold on
                    num_neigh_i=round(num_of_neighs(j));
                    if num_neigh_i<1
                        num_neigh_i=1;
                    end
                    plot(trj_2_red(j,1)*obj.scale,trj_2_red(j,2)*obj.scale,'o','Color',cmap(num_neigh_i,:))
                    xlabel('X (\mum)')
                    ylabel('Y (\mum)')
                    title(['Timestep ' num2str(dti)])
                    hold off
                end
                line([0 0],[-1000 1000],'Color','Black');
                line([-1000 1000],[0 0],'Color','Black');
                axis equal
                axis(axislimits)
                
            end
            
        end
        function output=find_reduced_trajectory_temp(obj,dt,noshow)
            %This function gives the tracks around the frame of reference
            %of particle each central particle
            %The input variable noshow indicates whether or not to plot the
            %reduced tracks. Give (1) for yes and (0) for no. Default is
            %no.
            %dt is the stepsize over which you want to measure the reduced
            %trajectory. Default is 1
            if nargin<3
                noshow=0;
            end
            if nargin<2
                dt=1;
                noshow=0;
            end
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            trace=obj.tr;
            nof=obj.NOF;
            if nof<3
                disp('There are fewer than 3 frames. That is not enough to determine the rotation angle.')
                output=NaN;
            else
                output={};
                %For each particle
                for p=1:obj.NOP
                    %Define three trajectories with a difference of 1
                    %timestep
                    trace_p=trace(find(trace(:,4)==p),1:2);
                    trj_0=trace_p(1+dt:dt:nof-dt,1:2);
                    trj_1=trace_p(1:dt:nof-2*dt,1:2);
                    trj_2=trace_p(1+2*dt:dt:nof,1:2);
                    
                    %Define the reduced trajectories
                    trj_0_red=zeros(length(trj_0(:,1)),2);
                    trj_1_red=trj_1(:,1:2)-trj_0(:,1:2);
                    trj_2_red=trj_2(:,1:2)-trj_0(:,1:2);
                    
                    %Rotate the reduced trajectories to align with the
                    %x-axis
                    trj_1_rot=[];
                    trj_2_rot=[];
                    for j=1:length(trj_0_red)
                        angle=pi-atan2(trj_1_red(j,2),trj_1_red(j,1));
                        R=[cos(angle)  -sin(angle) ; sin(angle)  cos(angle)];
                        trj_1_rot(j,1:2)=R*trj_1_red(j,1:2)';
                        trj_2_rot(j,1:2)=R*trj_2_red(j,1:2)';
                    end
                    if noshow==1
                        figure
                        hold on
                        box on
                        axis equal
                        xlabel('X (\mum)')
                        ylabel('Y (\mum)')
                        plot(trj_0_red(:,1)*obj.scale,trj_0_red(:,2)*obj.scale,'ko')
                        plot(trj_1_rot(:,1)*obj.scale,trj_1_rot(:,2)*obj.scale,'r-')
                        plot(trj_2_rot(:,1)*obj.scale,trj_2_rot(:,2)*obj.scale,'b-')
                        hold off
                    end
                    tr_red_temp=[trj_0_red(:,1:2) ones(length(trj_0_red(:,1)),1)*0 ; trj_1_rot(:,1:2) ones(length(trj_1_rot(:,1)),1)*1 ; trj_2_rot(:,1:2) ones(length(trj_2_rot(:,1)),1)*2];
                    output{p}=tr_red_temp;
                    
                end
            end
        end
        
        function obj=correlate_motion(obj,dt_low,dt_step,dt_high,ID)
            %Takes the distplacement histograms of a number of timesteps dt
            %And returns a graph of the uncorrelated motion, correlated
            %motion and angle in the motion with time. From those slopes,
            %the diffusion coefficient and instantaneous velocity can be
            %obtained with high precision
            
            if nargin==1
                dt_low=1;
                dt_step=1;
                dt_high=10;
            end
            if nargin==2
                dt_step=1;
                dt_low=round(dt_low);
                if dt_low<1
                    dt_low=1;
                end
                dt_high=dt_low+10;
            end
            if nargin==3
                dt_step=round(dt_step);
                if dt_step<1
                    dt_step=1;
                end
                dt_low=round(dt_low);
                if dt_low<1
                    dt_low=1;
                end
                dt_high=dt_low+10*dt_step;
            end
            if nargin>3
                dt_step=round(dt_step);
                if dt_step<1
                    dt_step=1;
                end
                dt_low=round(dt_low);
                if dt_low<1
                    dt_low=1;
                end
                dt_high=round(dt_high);
                if dt_high<dt_low+dt_step
                    dt_high=dt_low+dt_step;
                end
            end
            if nargin==5
                p_end=ID;
                p_begin=ID;
            else
                p_end=obj.NOP;
                p_begin=1;
            end
            
            nsteps=round((dt_high-dt_low)/dt_step);
            
            
            v0=NaN(nsteps+1,2);
            alpha0=NaN(nsteps+1,1);
            D0=NaN(nsteps+1,1);
            v0(1,:)=0;
            alpha0(1,1)=0;
            D0(1,1)=0;
            count=1;
            output={};
            for dt=dt_low:dt_step:dt_high
                count=count+1;
                output{count}=find_reduced_trajectory_temp(obj,dt);
                
                trj_0_red=[];
                trj_1_red=[];
                trj_2_red=[];
                for p=p_begin:p_end
                    outputi=output{count};
                    tr_red_p=outputi{p};
                    %find trajectories of particle p
                    trj_0_red_p=tr_red_p(find(tr_red_p(:,3)==0),1:2);
                    trj_1_red_p=tr_red_p(find(tr_red_p(:,3)==1),1:2);
                    trj_2_red_p=tr_red_p(find(tr_red_p(:,3)==2),1:3);
                    
                    %append into trajectories of all particle
                    trj_0_red=[trj_0_red ; trj_0_red_p];
                    trj_1_red=[trj_1_red ; trj_1_red_p];
                    trj_2_red=[trj_2_red ; trj_2_red_p];
                end
                v0x=nanmean(trj_2_red(:,1));
                v0y=nanmean(trj_2_red(:,2));
                v0(count,1)=sqrt(v0x^2+v0y^2);
                v0(count,2)=dt;
                alpha0(count,1)=atan2(v0y,v0x);
                
                trj_2_red_cent=[trj_2_red(:,1)-v0x trj_2_red(:,2)-v0y trj_2_red(:,3)];
                dist2=sqrt(trj_2_red_cent(:,1).^2+trj_2_red_cent(:,2).^2);
                %histdist=histogram(dist2);
                %figure
                %plot(histdist)
                %info=histfit(dist2);
                %D0(dt,1)=info.sigma;
                D0(count,1)=nanstd(dist2);
            end
            
            %Plot the bare velocity, diffusion and angle in time
            figure
            hold on
            box on
            xlabel(' time (frame) ','FontSize',14)
            ylabel(' v_0 (pixel)  ','FontSize',14)
            plot(v0(:,2),v0(:,1),'ko','DisplayName','Persistent speed');
            plot(v0(:,2),D0(:,1),'ro','DisplayName','Translational diffusion');
            %plot(v0(:,2),alpha0(:,1),'bo');
            hold off
            
            %Fit the correlated displacement to find the velocity
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',0,...
                'Upper',Inf,...
                'StartPoint',1);
            ft = fittype('v*x','options',fo);
            
            f1=fit(v0(:,2),v0(:,1),ft);
            
            figure
            hold on
            box on
            title('Correlated / swimming motion')
            xlabel(' time (frame) ','FontSize',14)
            ylabel(' v_0 (pixel)  ','FontSize',14)
            plot(v0(:,2),v0(:,1),'ko');
            plot(linspace(0.8*min(v0(:,2)),1.2*max(v0(:,2)),100),f1.v*linspace(0.8*min(v0(:,2)),1.2*max(v0(:,2)),100),'k-');
            hold off
            
            %Fit the uncorrelated displacement to find the diffusion
            %coefficient
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Lower',0,...
                'Upper',Inf,...
                'StartPoint',1);
            ft = fittype('4*Diff*x','options',fo);
            
            f2=fit(v0(:,2),D0(:,1).^2,ft);
            
            figure
            hold on
            box on
            title('Diffusive motion')
            xlabel(' time (frame) ','FontSize',14)
            ylabel(' MSD (pixel)  ','FontSize',14)
            plot(v0(:,2),(D0(:,1)).^2,'ko');
            plot(linspace(0.8*min(v0(:,2)),1.2*max(v0(:,2)),100),4*f2.Diff*linspace(0.8*min(v0(:,2)),1.2*max(v0(:,2)),100),'k-');
            hold off
            
            obj.motion=[f1.v f2.Diff];
            
            
        end
        function obj=crop_in_time(obj,tbegin,tend)
            %Returns the object with only data in between tbegin and tend
            %The only items it crops are obj.pos and obj.tr;
            
            if isempty(obj.pos)
                obj=find_pos(obj);
            end
            
            if tbegin<1
                disp('Begin time cannot be smaller than 1.')
                tbegin=1;
            end
            
            if ~isempty(obj.trparam)
                if obj.trparam{1,1}.good>(tend-tbegin)
                    %If the new size of the object is smaller than the required
                    %length for an accepted track, shorten that required
                    %length.
                    obj.trparam{1,1}.good=(tend-tbegin);
                end
            end
            
            nbegin=min(find(obj.pos(:,3)==tbegin));
            nend=max(find(obj.pos(:,3)==tend));
            obj.pos=obj.pos(nbegin:nend,:);
            if ~isempty(obj.tr)
                obj=track_obj(obj,0);
            end
            
        end
        function obj=calc_speed_vs_time(obj,plotcheck,particlelist)
            %Returns a graph with the instantaneous speed for each particle
            %in time
            %particlelists needs to be a 1xn vector listing the particle
            %IDS
            %to plot in the speeds vs time graph
            
            dt=8;
            
            if nargin==1
                particlelist=linspace(1,obj.NOP,obj.NOP);
                plotcheck=0;
            end
            
            if nargin==2
                particlelist=linspace(1,obj.NOP,obj.NOP);
            end
            
            %Assert that [tr] matrix exists
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            %get particle traces
            trace=obj.tr;
            %create empty matrix to save particle velocities
            spe=NaN(obj.NOF-dt,obj.NOP);
            %for each particle, calculate the velocity vector
            for i=1:obj.NOP
                %stagger the position matrices to get the velocity vector
                %in a single operation
                tr_i=trace(find(trace(:,4)==i),1:3);
                tr_x1=tr_i(1+dt:length(tr_i),1);
                tr_x2=tr_i(1:length(tr_i)-dt,1);
                tr_y1=tr_i(1+dt:length(tr_i),2);
                tr_y2=tr_i(1:length(tr_i)-dt,2);
                spe_tempx=tr_x1-tr_x2;
                spe_tempy=tr_y1-tr_y2;
                spe(:,i)=1/dt*sqrt(spe_tempx.^2+spe_tempy.^2);
            end
            obj.speeds=spe;
            
            if plotcheck==1
                %Get color scheme jet for length of matrix size spe.
                cmap=jet(length(particlelist));
                %plot the instantaneous speed in a figure
                figure
                hold on
                box on
                title('Instantaneous speeds')
                xlabel(' Time (s) ','FontSize',14)
                ylabel(' Speed (\mum/s)  ','FontSize',14)
                for q=1:length(particlelist)
                    i=particlelist(q);
                    plot(tr_i(1+dt:length(tr_i),3)/obj.framerate,spe(:,i)*obj.scale*obj.framerate,'Color',cmap(q,:),'DisplayName',['Particle ' num2str(i)])
                end
                %lgd=legend('show');
                hold off
                
                
                %get mean speed
                spe_mean=nanmean(spe,2);
                
                figure
                hold on
                box on
                title('Average speed all particles')
                xlabel(' Time (s) ','FontSize',14)
                ylabel(' Speed (\mum/s)  ','FontSize',14)
                plot(tr_i(1+dt:length(tr_i),3)/obj.framerate,spe_mean(:,1)*obj.scale*obj.framerate,'ko','DisplayName','Average');
                hold off
            end
            average_speed=nanmean(spe_mean(:,1));
            std_speed=nanstd(spe_mean(:,1));
            disp(['The average speed is ' num2str(average_speed) '.']);
            disp(['The standard deviation in the speed is ' num2str(std_speed) '.']);
        end
        function obj=calc_speed_vs_time_interval(obj,dt,particlelist)
            %Returns a graph with the instantaneous speed for each particle
            %in time
            %particlelists needs to be a 1xn vector listing the particle
            %IDS
            %to plot in the speeds vs time graph
            
            if nargin==1
                dt=1;
                particlelist=linspace(1,obj.NOP,obj.NOP);
            elseif nargin==2
                particlelist=linspace(1,obj.NOP,obj.NOP);
            end
            
            %Assert that [tr] matrix exists
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            %get particle traces
            trace=obj.tr;
            %create empty matrix to save particle velocities
            spe=NaN(obj.NOF-dt,obj.NOP);
            %for each particle, calculate the velocity vector
            for i=1:obj.NOP
                %stagger the position matrices to get the velocity vector
                %in a single operation
                tr_i=trace(find(trace(:,4)==i),1:3);
                tr_x1=tr_i(1+dt:length(tr_i),1);
                tr_x2=tr_i(1:length(tr_i)-dt,1);
                tr_y1=tr_i(1+dt:length(tr_i),2);
                tr_y2=tr_i(1:length(tr_i)-dt,2);
                spe_tempx=tr_x1-tr_x2;
                spe_tempy=tr_y1-tr_y2;
                spe(:,i)=sqrt(spe_tempx.^2+spe_tempy.^2)./dt;
            end
            obj.speeds=spe;
            
            %Get color scheme jet for length of matrix size spe.
            cmap=jet(length(particlelist));
            %plot the instantaneous speed in a figure
            figure
            hold on
            box on
            title('Instantaneous speeds')
            xlabel(' Time (s) ','FontSize',14)
            ylabel(' Speed (\mum/s)  ','FontSize',14)
            for q=1:length(particlelist)
                i=particlelist(q);
                plot(linspace(1,length(spe(:,1)),length(spe(:,1)))./obj.framerate,spe(:,i)*obj.scale*obj.framerate,'Color',cmap(q,:),'DisplayName',['Particle ' num2str(i)])
            end
            %lgd=legend('show');
            hold off
            
            %get mean speed
            spe_mean=nanmean(spe,2);
            
            figure
            hold on
            box on
            title('Average speed all particles')
            xlabel(' Time (s) ','FontSize',14)
            ylabel(' Speed (\mum/s)  ','FontSize',14)
            plot(linspace(1,length(spe_mean(:,1)),length(spe_mean(:,1)))./obj.framerate,spe_mean(:,1)*obj.scale*obj.framerate,'ko','DisplayName','Average');
            hold off
            
        end
        function obj=calc_speed_vs_time_neighbours(obj,ID)
            %Returns a graph with the instantaneous speed for each particle
            %in time
            %particlelists needs to be a 1xn vector listing the particle
            %IDS
            %to plot in the speeds vs time graph
            %speeds is returned as a 2D matrix containing time on the
            %verticle axis and particle ID on the horizontal axis
            
            
            %Assert that [tr] matrix exists
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            if isempty(obj.num_neighs)
                obj=find_neighbours(obj);
            end
            
            particlelist=[ID obj.neighbours{ID}];
            
            %get particle traces
            trace=obj.tr;
            %create empty matrix to save particle velocities
            spe=NaN(obj.NOF-1,obj.NOP);
            %for each particle, calculate the velocity vector
            for i=1:obj.NOP
                %stagger the position matrices to get the velocity vector
                %in a single operation
                tr_i=trace(find(trace(:,4)==i),1:3);
                tr_x1=tr_i(2:length(tr_i),1);
                tr_x2=tr_i(1:length(tr_i)-1,1);
                tr_y1=tr_i(2:length(tr_i),2);
                tr_y2=tr_i(1:length(tr_i)-1,2);
                spe_tempx=tr_x1-tr_x2;
                spe_tempy=tr_y1-tr_y2;
                spe(:,i)=sqrt(spe_tempx.^2+spe_tempy.^2);
            end
            obj.speeds=spe;
            
            %Get color scheme jet for length of matrix size spe.
            cmap=winter(length(particlelist));
            %plot the instantaneous speed in a figure
            figure
            hold on
            box on
            title('Instantaneous speeds')
            xlabel(' Time (s) ','FontSize',14)
            ylabel(' Speed (\mum/s)  ','FontSize',14)
            for q=1:length(particlelist)
                i=particlelist(q);
                plot(tr_i(2:length(tr_i),3),spe(:,i),'Color',cmap(q,:),'DisplayName',['Particle ' num2str(i)])
            end
            %lgd=legend('show');
            hold off
            
        end
        function obj=find_neighbours_1(obj,cutoffvalue)
            %creates matrix num_neighs(indicating how many neighbours each particle has)
            %creates matrix neighbours (indicating which particles are the
            %neighbours)
            %neighbours are defined using a cutoff value
            %Assumes that the number of neighbours is a fixed value for one
            %object and does not change from frame to frame
            %uses average distance between particles to find neighbours
            
            
            if nargin<2
                if ~isempty(obj.cutoff)
                    cutoffvalue=obj.cutoff;
                else
                    cutoffvalue=20;
                end
            end
            
            %if dist_av does not exist, make dist_av matrix first
            if isempty(obj.dist_av)
                if isempty(obj.dist)
                    obj=make_av_dist_matrix(obj);
                    av_dist=obj.dist_av;
                else
                    distance=obj.dist;
                    av_dist=nanmean(distance,3);
                end
            else
                av_dist=obj.dist_av;
            end
            
            nops=length(av_dist(:,1));
            num_neighs_temp=NaN(nops,1);
            obj.neighbours={};
            for i=1:nops
                neighbourstest=find(av_dist(i,:)<cutoffvalue & av_dist(i,:)>0);
                neighbours_temp=neighbourstest(find(neighbourstest~=i));
                obj.neighbours{i}=neighbours_temp;
                num_neighs_temp(i,1)=length(neighbourstest)-1;
            end
            obj.num_neighs=num_neighs_temp;
            obj.cutoff=cutoffvalue;
            disp(['The largest number of neighbours found was ' num2str(max(num_neighs_temp)) '.'])
        end
        function obj=find_neighbours_time(obj,cutoffvalue)
            %Finds number of neighbours for each particle in each frame
            %and saves in a 2D matrix.
            % times goes in width, ID in length.
            
            %Make sure a proper value of cutoff is used
            if nargin<2
                if ~isempty(obj.cutoff_time)
                    cutoffvalue=obj.cutoff_time;
                else
                    cutoffvalue=20;
                end
            end
            
            %if dist does not exist, make dist matrix first
            if isempty(obj.dist)
                obj=make_dist_matrix(obj);
            end
            
            %Get the distance matrix, number of particles and number of
            %frames
            distance=obj.dist;
            nops=obj.NOP;
            nofs=obj.NOF;
            
            %Create an 3D empty matrix for the number of neighbours
            neighbour=NaN(nops,nofs);
            
            %loop over all frames
            for t=1:nofs
                %loop over all particle IDS
                for i=1:nops
                    %calculate the number of neighbours for each i&t
                    neighbourstest=find(distance(i,:,t)<cutoffvalue);
                    %save number of neighbours to empy matrix
                    if ~isempty(neighbourstest)
                        neighbour(i,t)=length(neighbourstest)-1;
                    else
                        neighbour(i,t)=NaN;
                    end
                end
            end
            
            %save the obtained values to the object
            obj.numneighs_time=neighbour;
            obj.cutoff_time=cutoffvalue;
            disp(['The largest number of neighbours found was ' num2str(max(max(neighbour))) '.'])
            
            %get a colormap going through all particle IDs
            cmap=jet(nops);
            %plot the number of neighbours in time
            figure
            hold on
            box on
            title('Number of neighbours')
            xlabel(' Time (s) ','FontSize',14)
            ylabel(' Number of neighbours ','FontSize',14)
            for i=1:nops
                plot(linspace(1,nofs,nofs),neighbour(i,:),'o','Color',cmap(i,:),'DisplayName',['Particle ' num2str(i)]);
            end
            lgd=legend('show');
            hold off
        end
        function correlate_speed_numneighs(obj,cutoffvalue,particlelist)
            %plots the instantaneous speed of a particle versus its
            %instantaneous number of neighbours
            %particle needs to be a 1xn list of IDs
            
            %Define default cutoffvalue if not provided by user
            if nargin<2
                cutoffvalue=20;
            end
            if nargin<3
                if isempty(obj.tr)
                    obj=track_obj(obj);
                end
                particlelist=linspace(1,obj.NOP,obj.NOP);
            end
            
            %Assert that number of neighbours and speeds were previously
            %calculated and saved in the object
            if isempty(obj.numneighs_time)
                obj=find_neighbours_time(obj,cutoffvalue);
            end
            if isempty(obj.speeds)
                obj=calc_speed_vs_time(obj);
            end
            
            %define number of particles and a colourmap of that size
            nops=obj.NOP;
            spe=obj.speeds;
            numneighs=obj.numneighs_time;
            cmap=jet(nops);
            
            %create the matrix that correlates number of neighbours and
            %instantaneous speed
            corr=NaN(obj.NOF-1,2,nops);
            for q=1:length(particlelist)
                p=particlelist(q);
                for t=1:obj.NOF-1
                    corr(t,1,p)=numneighs(p,t);
                    corr(t,2,p)=spe(t,p);
                end
            end
            
            %average numneighs
            maxneigh=max(max(numneighs));
            y=NaN(1,maxneigh);
            for i=1:maxneigh
                y_temp=NaN(1,length(corr(1,1,:)));
                for p=1:length(corr(1,1,:))
                    corr_temp=corr(:,:,p);
                    y_temp(1,p)=nanmean(corr_temp(find(corr_temp(:,1)==i),2));
                end
                y(1,i)=nanmean(y_temp);
            end
            
            %plot figure
            figure
            hold on
            box on
            title('Correlation speed and aggregation number')
            xlabel(' Number of neighbours ','FontSize',14)
            ylabel(' Instantaneous speed (\mum/s ','FontSize',14)
            for i=1:nops
                plot(corr(:,1,i),corr(:,2,i),'-','Color',cmap(i,:),'DisplayName',['Particle ' num2str(i)]);
            end
            hold off
            
            figure
            hold on
            box on
            title('Average speed')
            xlabel(' Number of neighbours ','FontSize',14)
            ylabel(' Instantaneous speed (\mum/s ','FontSize',14)
            plot(linspace(1,maxneigh,maxneigh),y,'-');
            hold off
            
        end
        function output=corr_xy_motion(obj,selection,imshow)
            %for individual particles (nneighs==0)
            %plot the xx, yy, xy and yx correlation functions
            %selection can be any subset of all particles organized as [ID1 ID2 ID3...]
            %default selection is all particles
            %if selection=0, all particles are shown
            %imshow indicates whether (1) or not (0) to show graphs of
            %individual particles. Default is yes (1);
            
            if nargin<2
                selection=1:obj.NOP;
            end
            if selection==0
                selection=1:obj.NOP;
            end
            if nargin<3
                imshow=1;
            end
            
            %Check if number of neighbours is calculated
            if isempty(obj.num_neighs)
                obj=find_neighbours(obj);
            end
            
            %temporarily save number of neighbours
            neigh_temp=obj.num_neighs;
            
            if min(neigh_temp)>0
                disp('All particles have 1 or more neighbours and no individual correlation can be calculated.')
                output=[];
            else
                %For all particles that have no neighbours
                output=cell(1,3);
                count=1;
                for p=1:obj.NOP
                    corrs=cell(4,1);
                    if neigh_temp(p)==0 && any(selection==p)
                        %find trajectory of particle i
                        tr_temp=obj.tr(find(obj.tr(:,4)==p),1:3);
                        %define x and y coordinates in time
                        x=[tr_temp(:,3) tr_temp(:,1)];
                        y=[tr_temp(:,3) tr_temp(:,2)];
                        %correlate auto and cross
                        corrs{1,1}=correlate_functions(x,x);
                        corrs{2,1}=correlate_functions(x,y);
                        corrs{3,1}=correlate_functions(y,x);
                        corrs{4,1}=correlate_functions(y,y);
                        
                        %plot correlation functions
                        name2={'dx.dx (\mum^2/s)' 'dx.dy (\mum^2/s)' 'dy.dx (\mum^2/s)' 'dy.dy (\mum^2/s)'};
                        figure(30)
                        title(['Particle ' num2str(p)])
                        Dtens=NaN(2,2);
                        for i=1:4
                            axi=subplot(2,2,i);
                            %find fit parameter
                            x=corrs{i}(:,1)/obj.framerate;
                            y=corrs{i}(:,2)*obj.scale^2;
                            b=x\y;
                            axis(axi,[0 1.1*max(corrs{1}(:,1)/obj.framerate) -1.2*max(corrs{1}(:,2)*obj.scale^2) 1.2*max(corrs{1}(:,2)*obj.scale^2)]);
                            %errorbar(corrs{i}(:,1)/obj.framerate,corrs{i}(:,2)*obj.scale^2/obj.framerate,corrs{i}(:,3)*obj.scale^2/obj.framerate,'ko')
                            hold on
                            box on
                            plot(axi,x,y,'k.')
                            plot(axi,x,b*x,'r-')
                            hold off
                            %title(['Particle ' num2str(p) ', ' name{i};])
                            xlabel(axi,'dt (s)')
                            ylabel(axi,name2{i});
                            if i<3
                                Dtens(1,i)=b;
                            else
                                Dtens(2,i-2)=b;
                            end
                        end
                        if imshow~=1
                            close 30
                        end
                        %save findings in output
                        output{count,1}=p;
                        output{count,2}=corrs;
                        output{count,3}=Dtens;
                        count=count+1;
                    end
                    
                end
            end
            
            %Average the correlation and plot average
            allxx=NaN(length(output{1,2}{1,1}(:,1)),length(output(:,1)));
            allxy=NaN(length(output{1,2}{1,1}(:,1)),length(output(:,1)));
            allyx=NaN(length(output{1,2}{1,1}(:,1)),length(output(:,1)));
            allyy=NaN(length(output{1,2}{1,1}(:,1)),length(output(:,1)));
            for q=1:length(output(:,1))
                allxx(:,q)=output{q,2}{1,1}(:,2);
                allxy(:,q)=output{q,2}{2,1}(:,2);
                allyx(:,q)=output{q,2}{3,1}(:,2);
                allyy(:,q)=output{q,2}{4,1}(:,2);
            end
            meanall=cell(4,1);
            meanall{1,1}=nanmean(allxx,2);
            meanall{2,1}=nanmean(allxy,2);
            meanall{3,1}=nanmean(allyx,2);
            meanall{4,1}=nanmean(allyy,2);
            x2=linspace(1,obj.NOF/100,obj.NOF/100)';
            
            
            %plot average correlation function
            name2={'dx.dx (\mum^2)' 'dx.dy (\mum^2)' 'dy.dx (\mum^2)' 'dy.dy (\mum^2)'};
            f = figure;
            p = uipanel('Parent',f,'BorderType','none','BackgroundColor','white');
            p.Title = 'Average';
            p.TitlePosition = 'centertop';
            p.FontSize = 12;
            p.FontWeight = 'bold';
            Dav=NaN(2,2);
            for i=1:4
                axi=subplot(2,2,i,'Parent',p);
                %find fit parameter
                x=x2(:,1)/obj.framerate;
                y=meanall{i,1}*obj.scale^2;
                b=x\y;
                axis(axi,[0 1.1*max(x(:,1)) -1.2*max(meanall{1,1}(:,1)*obj.scale^2) 1.2*max(meanall{1,1}(:,1)*obj.scale^2)])
                %errorbar(corrs{i}(:,1)/obj.framerate,corrs{i}(:,2)*obj.scale^2/obj.framerate,corrs{i}(:,3)*obj.scale^2/obj.framerate,'ko')
                hold on
                box on
                plot(axi,x,y,'k.')
                plot(axi,x,b*x,'r-')
                hold off
                %title(['Particle ' num2str(p) ', ' name{i};])
                xlabel(axi,'dt (s)')
                ylabel(axi,name2{i});
                if i<3
                    Dav(1,i)=b;
                else
                    Dav(2,i-2)=b;
                end
            end %end forloop image
            output{count,1}='Average';
            output{count,2}=meanall;
            output{count,3}=Dav;
        end %end function
        function calc_trace_distance_histogram(obj)
            %Takes the trajectory of a particle and calculates the distance
            %of each point in the trajectory with each other point for
            %various effective framerates (minimal timesteps).
            %produces a histogram of this information
            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            distance=cell(1,1);
            for p=1:obj.NOP
                trp=obj.tr(find(obj.tr(:,4)==p),1:2);
                dista=[];
                for dt=5:length(trp-1)
                    x1=trp(1:end-dt,1);
                    x2=trp(1+dt:end,1);
                    y1=trp(1:end-dt,2);
                    y2=trp(1+dt:end,2);
                    distp=sqrt((x2-x1).^2+(y2-y1).^2)*obj.scale;
                    dista=[dista ;distp];
                end
                distance{p}=dista;
            end
            
            figure
            count1=floor(sqrt(obj.NOP))+1;
            count2=floor(sqrt(obj.NOP))+1;
            for p=1:obj.NOP
                if sum(~isnan(distance{p}))>1
                    subplot(count1,count2,p)
                    hold on
                    box on
                    histogram(distance{p},'FaceColor','b','FaceAlpha',1)
                    xlabel('Stepize (\mum)');
                    ylabel('N');
                    axis([0 1.2*max(distance{p}) -1 1]);
                    axis 'auto y'
                    hold off
                end
            end
            
        end
        function obj=smooth_trace(obj,plotcheck)
            %smoothes the trace of object using a local regression using
            %weighted linear least squares and a 1st degree polynomial model
            
            if nargin==1
                plotcheck=0;
            end
            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            for p=1:obj.NOP
                trp=obj.tr(find(obj.tr(:,4)==p),:);
                trp_smooth=trp;
                trp_smooth(:,1)=smooth(trp(:,3),trp(:,1),19);
                trp_smooth(:,2)=smooth(trp(:,3),trp(:,2),19);
                
                if plotcheck==1
                    figure
                    hold on
                    plot(trp(:,3),trp(:,1),'r-')
                    
                    plot(trp_smooth(:,3),trp_smooth(:,1),'b-')
                    hold off
                    
                    figure
                    hold on
                    plot(trp(:,3),trp(:,2),'r-')
                    
                    plot(trp_smooth(:,3),trp_smooth(:,2),'b-')
                    hold off
                end
                obj.tr(find(obj.tr(:,4)==p),:)=trp_smooth;
            end
            
            
        end %end function
        function obj=manual_tracking(obj,n,dt)
            %finds positions by manual particle finding
            %n is number of particles you wish to find in each frame
            
            if nargin<2
                n=1;
            end
            
            fname=[obj.foldername obj.filename '.tif'];
            
            
            
            info = imfinfo(fname);
            num_images = numel(info);
            pos_temp=NaN(n*num_images,3);
            
            for p=1:dt:num_images
                picture=imread(fname,p);
                figure(p)
                colormap('gray'),imagesc(picture)
                truesize(p,[700 700])
                [x,y]=ginput(n);
                close(p);
                if n==1
                    pos_temp(p,1:3)=[x y p];
                else
                    pos_temp(n*(p-1)+1:n*p,1:3)=[x y p*ones(n,1)];
                end
                
            end
            
            obj.pos=pos_temp;
            
        end
        function obj=manualy_fix_tracking(obj,beginframe,endframe)
            %finds positions by manual particle finding
            %m is number of frames you want to be doing this for
            
            if isempty(obj.pos)
                obj=find_pos(obj);
            end
            
            fname=[obj.foldername obj.filename '.tif'];
            
            info = imfinfo(fname);
            if nargin<2
                beginframe=1;
                num_images = numel(info);
            elseif nargin<3
                beginframe=1;
                num_images = endframe;
            else
                num_images=endframe;
            end
            
            
            pos_old=obj.pos;
            if beginframe~=1
                pos_new=pos_old(find(pos_old(:,3)<beginframe),1:3);
            else
                pos_new=[];
            end
            
            for p=beginframe:num_images
                picture=imread(fname,p);
                figure(p)
                pos_old_p=pos_old(find(pos_old(:,3)==p),1:3);
                hold on
                colormap('gray'),imagesc(picture)
                viscircles(pos_old_p(:,1:2),ones(length(pos_old_p(:,1)),1)*2)
                truesize(p,[700 700])
                hold off
                try n=get_numerical_input('How many particles need be added?');
                catch
                    try n=get_numerical_input('How many particles need be added?');
                    catch
                        n=0;
                    end
                end
                if n~=0
                    try [x,y]=ginput(n);
                    catch
                        x=NaN;
                        y=NaN;
                    end
                    pos_temp=[x y p*ones(length(y(:,1)),1)];
                    pos_new_p=[pos_old_p ; pos_temp];
                else
                    pos_new_p=pos_old_p;
                end
                
                pos_new=[pos_new ; pos_new_p];
                try close(p);
                catch
                end
            end
            
            if num_images~=numel(info)
                pos_after=pos_old(find(pos_old(:,3)>num_images),1:3);
            end
            
            pos_new=[pos_new ; pos_after];
            obj.pos=pos_new;
            
        end
        function obj_new=split_trace_new(obj,p)
            %takes one particle and places it in a new object
            
            params=get_parameters(obj);
            fname=[obj.foldername obj.filename '.tif'];
            obj_new=particles(fname,params);
            obj_new.scale=obj.scale;
            obj_new.framerate=obj.framerate;
            postemp=obj.tr(find(obj.tr(:,4)==p),1:3);
            trtemp=[postemp ones(length(postemp(:,1)),1)];
            obj_new=add_trace(obj_new,trtemp);
            obj_new=add_pos(obj_new,postemp);
            
        end
        function obj=fix_pos_nans(obj)
            %ONLY WORKS ON OBJECT WITH ONE PARTICLE CREATED WITH
            %split_trace_new;
            %finds positions by manual particle finding
            %m is number of frames you want to be doing this for
            
            if isempty(obj.pos)
                obj=find_pos(obj);
            end
            
            if obj.NOP~=1
                error('This function will not work for NOP>1');
            end
            
            fname=[obj.foldername obj.filename '.tif'];
            
            info = imfinfo(fname);
            num_images=numel(info);
            
            for p=1:num_images
                pos_p=obj.pos(p,:);
                if isnan(pos_p(1,1))
                    picture=imread(fname,p);
                    figure(p)
                    hold on
                    colormap('gray'),imagesc(picture)
                    plot(obj.tr(:,1),obj.tr(:,2),'k-');
                    truesize(p,[700 700])
                    hold off
                    try [x,y]=ginput(1);
                    catch
                        x=NaN;
                        y=NaN;
                    end
                    pos_new_p=[x y p];
                    obj.pos(p,1:3)=pos_new_p;
                    try close(p);
                    catch
                    end
                end
            end
            
        end
        function obj=calc_g_of_r(obj,binsize)
            %takes as an input the position of all particles in all frames
            %and calculates the g(r) given a certain binsize in each frame.
            %data is stored in 'gr' with format [r gr std r gr std...] as
            %an r x 2 NOP matrix
            %data is stored in pixels
            
            %IMPORTANT. This code only works for tracked movies. If you
            %have a single image, copy it twice (so that you have three and can
            %track it) then track it and run this code. Then simply use the
            %g(r) of one of the two (identical) images.
            %the reason for this flaw it that the particles need to be
            %labeld for this code to work
            
            %if no binsizes is provided, defaults is 1 (pixel)
            if nargin==1
                binsize=1;
            end
            
            if isempty(obj.dist)
                obj=make_dist_matrix(obj);
            end
            
            maxval=length(obj.im)/3;
            nbins=ceil(maxval/binsize);
            
            imsize=size(obj.im);
            
            gr_temp=NaN(nbins,3*obj.NOF);
            %loop over all frames
            for t=1:obj.NOF
                dist_t=obj.dist(:,:,t);
                n_av=NaN(nbins,1);
                n_std=NaN(nbins,1);
                for r=1:nbins
                    limlast=(r-1)*binsize;
                    lim=r*binsize;
                    xlimmin=0+lim;
                    ylimmin=0+lim;
                    xlimmax=imsize(1)-lim;
                    ylimmax=imsize(2)-lim;
                    n_neighs=NaN(obj.NOP,1);
                    %scalefactor to correct for increasing area with
                    %increasing r
                    scalefactor=lim^2-limlast^2;
                    for p=1:obj.NOP
                        xpos=obj.tr((p-1)*obj.NOF+t,1);
                        ypos=obj.tr((p-1)*obj.NOF+t,1);
                        if xpos<xlimmax && ypos<ylimmax && xpos>xlimmin && ypos>ylimmin
                            n_neighs_old=length(find(dist_t(:,p)<=lim & dist_t(:,p)>limlast));
                            n_neighs(p,1)=n_neighs_old/scalefactor;
                        end
                    end
                    n_av(r,1)=nanmean(n_neighs(:,1));
                    n_std(r,1)=nanstd(n_neighs(:,1));
                end
                x=(1:nbins)*binsize-0.5*binsize;
                figure
                hold on
                box on
                errorbar(obj.scale*x',n_av,n_std,'ko-')
                hold off
                
                gr_temp(:,3*(t-1)+1)=x';
                gr_temp(:,3*(t-1)+2)=n_av/(nanmean(n_av(:,1)));
                gr_temp(:,3*(t-1)+3)=n_std/(nanmean(n_av(:,1)));
            end
            
            obj.gr=gr_temp;
        end
        function plot_g_of_r(obj,list,identity)
            %list is a list of times for which to plot the g_of_r
            %default is just the first image
            %identity indicates whether or not to plot errorbars
            %default (1) is yes
            
            if nargin<2
                list=1;
            end
            if nargin<3
                identity=1;
            end
            
            cmap=copper(length(list));
            count=1;
            
            if identity==1
                figure
                hold on
                box on
                for t=1:obj.NOF
                    if any(list)==t
                        errorbar(obj.scale*obj.gr(:,3*(t-1)+1), obj.gr(:,3*(t-1)+2),obj.gr(:,3*(t-1)+3),'o-','Color',cmap(count,:))
                        count=count+1;
                    end
                end
                axis([0 obj.scale*max(obj.gr(:,3*(t-1)+1)) 0 1+max(obj.gr(:,3*(t-1)+3))])
                xlabel('r (\mum)');
                ylabel('g(r)');
                hold off
            else
                figure
                hold on
                box on
                for t=1:obj.NOF
                    if any(list)==t
                        plot(obj.scale*obj.gr(:,3*(t-1)+1), obj.gr(:,3*(t-1)+2),'o-','Color',cmap(count,:))
                        count=count+1;
                    end
                end
                axis([0 obj.scale*max(obj.gr(:,3*(t-1)+1)) 0 3])
                xlabel('r (\mum)');
                ylabel('g(r)');
                hold off
            end
        end
    end %end methods
end % end classdef