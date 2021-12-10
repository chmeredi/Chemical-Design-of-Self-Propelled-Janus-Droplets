classdef dropletpair2<droplets
    %Dropletpairs is a subclass of droplets designed to investigate
    %interactions between two droplets
    %Droplets class is a subclass of particles
    %
    %
    %Workflow:
    %   - Open tiffstack with Droplets, find positions and track
    %   - run objnew=split_to_pair2(obj). Click CHASEE FIRST, CHASER SECOND 
    %   - run obj=crop_in_time_trace(obj,begint,endt)
    %   - run obj=smooth_trace(obj)
    %   - run obj=analyze(obj) and plot relevant quantities
    %
    %
    %
    %Extra properties
    %   - distance      %interparticle distance vs time
    %   - u_rel         %relative speed vs distance : derivative of interparticle distance 
    %   - u_sum         %total speed of pair vs distance: derivative of centerpos
    %   - centerpos     %position of center of mass of pair
    %   - u1            %speed of particle 1 tangenential to distance vs distance
    %
    %Extra methods
    %
    %   Calculation methods
    %   - obj=dropletpair(fname,parameters,identity) initialization
    %   - obj=calc_distance(obj)
    %   - obj=calc_centerpos(obj)
    %   - obj=calc_u_rel(obj)
    %   - obj=calc_u_sum(obj)
    %   - obj=calc_ui(obj) 
    %   - obj=analyze(obj);
    %   - obj=shift_u_rel(obj,amount);
    %   - obj=shift_u_sum(obj,amount);
    %   - obj=shift_u1(obj,amount);
    %   - obj=shift_u2(obj,amount);
    %
    %   Display methods
    %   - plot_distance(obj)
    %   - plot_centerpos(obj)
    %   - plot_u_rel(obj)
    %   - plot_u_sum(obj)
    %   - plot_ui(obj)
    %Created on 19-07-18 by Pepijn Moerman
    %Last modified: -
    
    properties
        distance
        centerpos
        u_rel
        u_sum
        u1
        u2        
    end
    methods
        
        %-------------------CALCULATION METHODS-----------------------%
        
        function obj=dropletpair2(fname,parameters,identity)
            %Droplets needs to be an object of type particles
            %Identity indicates whether or not users is requested to find
            %positions and fill in scale.
            
            if nargin==1
                parameters=[];
            elseif nargin<1
                parameters=[];
                fname=[];
            end
            if nargin<3
                identity='y';
            end
            obj=obj@droplets(fname,parameters,identity);
        end
        function obj=calc_distance(obj)
            %Finds the distance between the particles in the object as
            %function of time
            
            %check if input exists
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            %get traces of particle 1 and 2
            tr1=obj.tr(1:obj.NOF,:);
            tr2=obj.tr(obj.NOF+1:2*obj.NOF,:);
            
            %calculate the distance between them
            r(:,1)=sqrt((tr1(:,1)-tr2(:,1)).^2+(tr1(:,2)-tr2(:,2)).^2);
            %r(:,1)=smooth(tr1(:,3),r(:,1));
            
            %save distance and time
            obj.distance=[tr1(:,3) r];
        end
        function obj=calc_centerpos(obj)
            %Finds the position of the center of mass of the two droplets
            %as a function of time
            %centerpos has form [t x_cm y_cm]
            
            %check if input exists
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            %get traces of particle 1 and 2
            tr1=obj.tr(1:obj.NOF,:);
            tr2=obj.tr(obj.NOF+1:2*obj.NOF);
            
            %calculate the center of mass
            xm=(tr1(:,1)+tr2(:,1))/2;
            ym=(tr1(:,2)+tr2(:,2))/2;
            
            %save position of center of mass vs time
            obj.centerpos=[tr1(:,3) xm ym];
            
        end
        function obj=calc_u_rel(obj)
            %Calculates the relative dispacement between the two droplets
            %which is a measure of the interaction strength
            
            %check if input exists
            if isempty(obj.distance)
                obj=calc_distance(obj);
            end
            
            %get the distance between particles
            r(:,1)=obj.distance(:,2);
            
            %smooth the displacement
            %r(:,1)=smooth(r(:,1));
                        
            %calculate change in distance between 2 frames
            dr(:,1)=r(2:end,:)-r(1:end-1,:);
            
            %save distance change between 2 frames as speed
            obj.u_rel=[r(1:end-1,1) dr(:,1)];
               
        end
        function obj=calc_u_sum(obj)
            %Calculates the displacement of the center of mass parallel to
            %the axis that connects the two droplets
            
            %check if input exists
            if isempty(obj.centerpos)
                obj=calc_centerpos(obj);
            end
            
            %get the center of mass
            cm(:,1:2)=obj.centerpos(:,2:3);
            
            %smooth center of mass
            %cm(:,1)=smooth(cm(:,1));
            %cm(:,2)=smooth(cm(:,2));
            
            %get the positions of particle 1 and 2
            tr1=obj.tr(1:obj.NOF,1:2);
            tr2=obj.tr(obj.NOF+1:2*obj.NOF,1:2);
            
            %smooth the x and y displacement
            %tr1(:,1)=smooth(tr1(:,1));
            %tr1(:,2)=smooth(tr1(:,2));
            %tr2(:,1)=smooth(tr2(:,1));
            %tr2(:,2)=smooth(tr2(:,2));
                        
            %calculate the x and y displacement of the center of mass
            dx(:,1)=cm(2:end,1)-cm(1:end-1,1);
            dy(:,1)=cm(2:end,2)-cm(1:end-1,2);
            
            %calculate the total displacment of the center of mass
            drtot(:,1)=sqrt((dx(:,1)).^2+(dy(:,1)).^2);
            
            %calculate the angle of the center of mass displacement with the labframe
            cmangle(:,1)=atan2(dy(:,1),dx(:,1));
            
            %smooth angle
            %cmangle(:,1)=smooth(cmangle(:,1));
            
            %calculate the angle of the vector connecting the two droplets
            %with the labframe
            pangle(:,1)=atan2((tr2(1:end-1,2)-tr1(1:end-1,2)),(tr2(1:end-1,1)-tr1(1:end-1,1)));
            
            %smooth angle
            %pangle(:,1)=smooth(pangle(:,1));
            
            %calculate the displacement parallel to the axis connecting the
            %droplets by taking the inproduct            
            drpar=drtot.*abs(cos(cmangle-pangle));
            
            %save the parallel displacement
            obj.u_sum=[obj.distance(1:end-1,2) drpar];
            
        end
        function obj=calc_ui(obj)
            %calculates the displacment of droplet 1 parallel to the axis connecting
            %droplet 1 and droplet 2
            
            %check if input exists
            if isempty(obj.centerpos)
                obj=calc_centerpos(obj);
            end
            
            %get the position of particle 1 and 2 and the distance
            tr1(:,1:3)=obj.tr(1:obj.NOF,1:3);
            tr2(:,1:3)=obj.tr(obj.NOF+1:2*obj.NOF,1:3);
            dist(:,1)=obj.distance(:,2);
            
            %smooth traces
            %tr1(:,1)=smooth(tr1(:,3),tr1(:,1));
            %tr1(:,2)=smooth(tr1(:,3),tr1(:,2));
            %tr2(:,1)=smooth(tr2(:,3),tr2(:,1));
            %tr2(:,2)=smooth(tr2(:,3),tr2(:,2));
            
            %calculate the x and y displacement of droplet 1 and 2
            dx1(:,1)=tr1(2:end,1)-tr1(1:end-1,1);
            dy1(:,1)=tr1(2:end,2)-tr1(1:end-1,2);
            dx2(:,1)=tr2(2:end,1)-tr2(1:end-1,1);
            dy2(:,1)=tr2(2:end,2)-tr2(1:end-1,2);
            
            %calculate the total displacement of droplet 1 and 2
            dpos1(:,1)=sqrt((dx1).^2+(dy1).^2);
            dpos2(:,1)=sqrt((dx2).^2+(dy2).^2);
            
            %calculate the angle of the droplet displacement with the labframe
            cmangle1(:,1)=atan2(dy1(:,1),dx1(:,1));
            cmangle2(:,1)=atan2(dy2(:,1),dx2(:,1));
            
            %smooth angle
            %cmangle1(:,1)=smooth(cmangle1(:,1));
            %cmangle2(:,1)=smooth(cmangle2(:,1));
            
            %calculate the angle of the vector connecting the two droplets
            %with the labframe
            pangle(:,1)=atan2((tr2(1:end-1,2)-tr1(1:end-1,2)),(tr2(1:end-1,1)-tr1(1:end-1,1)));
            
            %smooth angle
            %pangle(:,1)=smooth(pangle(:,1));
            
            %calculate the displacement parallel to the axis connecting the
            %droplets by taking the inproduct
%             if abs(mean(cmangle1)-mean(pangle))<pi/2
%                 drpar1(:,1)=-dpos1(:,1).*cos(cmangle1(:,1)-pangle(:,1));
%                 mean(cmangle1)
%                 mean(pangle)
%                 'old'
%                 abs(mean(cmangle1)-mean(pangle))
%             else
%                 drpar1(:,1)=dpos1(:,1).*cos(cmangle1(:,1)-pangle(:,1));
%                 mean(cmangle1)
%                 mean(pangle)
%                 'new'
%                 abs(mean(cmangle1)-mean(pangle))
%             end
            %drpar1(:,1)=dpos1(:,1).*cos(cmangle1(:,1)-pangle(:,1));
            drpar1(:,1)=-dpos1(:,1).*cos(cmangle1(:,1)-pangle(:,1));
            drpar2(:,1)=dpos2(:,1).*cos(cmangle2(:,1)-pangle(:,1));
            
            %save the parallel displacement vs distance
            obj.u1=[dist(1:end-1,1) drpar1];
            obj.u2=[dist(1:end-1,1) drpar2];
            
        end
        function obj=crop_in_time_trace(obj,tbegin,tend)
            %saves only part of the trajectory between tbegin and tend
            %assumes only two particles in object
            
            %make sure proper input is provided
            if nargin<3
                error('Please provide a begin and end time.')
            end
            
            %make sure input is sensible
            if tbegin<1
                tbegin=1;
            end
            if tend>obj.NOF
                tend=obj.NOF;
            end
            if tbegin>tend
                error('Please make sure the begin time is smaller than the end time.')
            end
            
            %get trace particle 1 and crop
            tr1=obj.tr(1:obj.NOF,:);
            tr1=tr1(tbegin:tend,:);
            tr1(:,3)=linspace(1,length(tr1(:,3)),length(tr1(:,3)))';
            
            %get trace particle 2 and crop
            tr2=obj.tr(obj.NOF+1:2*obj.NOF,:);
            tr2=tr2(tbegin:tend,:);
            tr2(:,3)=linspace(1,length(tr1(:,3)),length(tr1(:,3)))';
            
            %combine traces
            trcombi=[tr1; tr2];
            
            %save traces to object
            obj.tr=trcombi;
            obj.NOF=length(tr1(:,1));
        end
        function obj=analyze(obj,tbegin,tend)
            %Combines the calculation functions for short notation:
            %   - obj=smooth_trace(obj) -- From particle class
            %   - obj=calc_distance(obj)
            %   - obj=calc_centerpos(obj)
            %   - obj=calc_u_rel(obj)
            %   - obj=calc_u_sum(obj)
            %   - obj=calc_ui(obj)
            %also incorporates obj=crop_in_time_trace if a tbegin AND tend
            %are provided. Otherwise this part is skipped
            
            %check if enough input is provided to crop tracks in time
            if nargin==3
                obj=crop_in_time_trace(obj,tbegin,tend);
            end
            
            %perform other analysis functions
            obj=calc_distance(obj);
            obj=calc_centerpos(obj);
            obj=calc_u_rel(obj);
            obj=calc_u_sum(obj);
            obj=calc_ui(obj);
            
        end
        function obj=shift_u_rel(obj,amount)
            %If u_rel does not decay to zero at large separation, this
            %function corrects by shifting the whole function a chosen
            %amount
            %amount is added to u_rel
            %give amount in micrometer/second
            
            %check if input exists
            if isempty(obj.u_rel)
                obj=calc_u_rel(obj);
            end
            
            obj.u_rel(:,2)=obj.u_rel(:,2)+amount/obj.scale/obj.framerate;
            
            plot_u_rel(obj);
            
        end
        function obj=shift_u_sum(obj,amount)
            %If u_rel does not decay to zero at large separation, this
            %function corrects by shifting the whole function a chosen
            %amount
            %amount is added to u_rel
            %give amount in micrometer/second
            
            %check if input exists
            if isempty(obj.u_sum)
                obj=calc_u_sum(obj);
            end
            
            obj.u_sum(:,2)=obj.u_sum(:,2)+amount/obj.scale/obj.framerate;
            
            plot_u_sum(obj);
            
        end
        function obj=shift_u1(obj,amount)
            %If u_rel does not decay to zero at large separation, this
            %function corrects by shifting the whole function a chosen
            %amount
            %amount is added to u_rel
            %give amount in micrometer/second
            
            %check if input exists
            if isempty(obj.u1)
                obj=calc_ui(obj);
            end
            
            obj.u1(:,2)=obj.u1(:,2)+amount/obj.scale/obj.framerate;
            
            plot_ui(obj);
            
        end
        function obj=shift_u2(obj,amount)
            %If u_rel does not decay to zero at large separation, this
            %function corrects by shifting the whole function a chosen
            %amount
            %amount is added to u_rel
            %give amount in micrometer/second
            
            %check if input exists
            if isempty(obj.u2)
                obj=calc_ui(obj);
            end
            
            obj.u2(:,2)=obj.u2(:,2)+amount/obj.scale/obj.framerate;
            
            plot_ui(obj);
            
        end
        
        %-------------------DISPLAY METHODS-----------------------%
        
        function plot_distance(obj)
            %plots the distance between two droplets as function of time
            
            %check if input exists
            if isempty(obj.distance)
                obj=calc_distance(obj);
            end
            
            %get time and distance
            t=obj.distance(:,1);
            r=obj.distance(:,2);
            
            %plot figure and provide scale and framerate
            disp(['Uses scale = ' num2str(obj.scale) ' and framerate ' num2str(obj.framerate) '.'])
            figure
            hold on
            plot(t(:,1)/obj.framerate,r(:,1)*obj.scale,'ko')
            xlabel('t (s)');
            ylabel('r (\mum)');
            pbaspect([1 1 1])
            axis([0 obj.NOF/obj.framerate 0 1.1*max(r(:,1))*obj.scale])
            box on
            hold off
            
        end
        function plot_centerpos(obj)
            %plots the trajectory of the center of mass of both droplets
            
            %check if input exists
            if isempty(obj.centerpos)
                obj=calc_centerpos(obj);
            end
            
            %get time and distance
            x=obj.centerpos(:,2);
            y=obj.centerpos(:,3);
            
            %create a colormap
            n=length(obj.centerpos(:,1))-1;
            cmap=jet(n);
            
            %plot figure
            figure
            hold on
            for i=1:n
                plot(x(i:i+1)*obj.scale,y(i:i+1)*obj.scale,'Color',cmap(i,:),'LineWidth',2)
            end
            xlabel('x (\mum)');
            ylabel('y (\mum)');
            pbaspect([1 1 1])
            box on
            hold off
                            
        end
        function plot_u_rel(obj)
            %plots the relative displacement between two particles
            
            %check if input exists
            if isempty(obj.u_rel)
                obj=calc_u_rel(obj);
            end
            
            %get time and distance
            r=obj.u_rel(:,1);
            u=obj.u_rel(:,2);
            
            %plot figure and provide scale and framerate
            disp(['Uses scale = ' num2str(obj.scale) ' and framerate ' num2str(obj.framerate) '.'])
            figure
            hold on
            plot(r(:,1)*obj.scale,u(:,1)*obj.scale*obj.framerate,'ko')
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            pbaspect([1 1 1])
            %axis([0 1.1*max(r(:,1))*obj.scale 1.1*min(u(:,1))*obj.scale*obj.framerate 1.3*max(u(:,1))*obj.scale*obj.framerate])
            box on
            hold off
            
        end
        function plot_u_sum(obj)
            %plots the displacement of the center of mass paralell to the axis connecting the two particles
            
            %check if input exists
            if isempty(obj.u_rel)
                obj=calc_u_rel(obj);
            end
            
            %get time and distance
            r=obj.u_sum(:,1);
            u=obj.u_sum(:,2);
            
            %plot figure and provide scale and framerate
            disp(['Uses scale = ' num2str(obj.scale) ' and framerate ' num2str(obj.framerate) '.'])
            figure
            hold on
            plot(r(:,1)*obj.scale,u(:,1)*obj.scale*obj.framerate,'ko')
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            pbaspect([1 1 1])
            %axis([0 1.1*max(r(:,1))*obj.scale 1.1*min(u(:,1))*obj.scale*obj.framerate 1.1*max(u(:,1))*obj.scale*obj.framerate])
            box on
            hold off
        end
        function plot_ui(obj)
            %plots the displacement of both droplets parallel to the axis connecting the two particles
            
            %check if input exists
            if isempty(obj.u_rel)
                obj=calc_u_rel(obj);
            end
            
            %get time and distance
            r(:,1)=obj.u1(:,1);
            v1(:,1)=obj.u1(:,2);
            v2(:,1)=obj.u2(:,2);
            
            %provide scale and framerate
            disp(['Uses scale = ' num2str(obj.scale) ' and framerate ' num2str(obj.framerate) '.'])
            
            %plot speed droplet 1
            figure
            hold on
            title('Particle 1')
            plot(r(:,1)*obj.scale,v1(:,1)*obj.scale*obj.framerate,'ko')
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            pbaspect([1 1 1])
            axis([0.9*min(r(:,1)) 1.1*max(r(:,1))*obj.scale 1.1*min(v1(:,1))*obj.scale*obj.framerate 1.1*max(v1(:,1))*obj.scale*obj.framerate])
            box on
            hold off
            
            %plot speed droplet 2
            figure
            hold on
            title('Particle 2')
            plot(r(:,1)*obj.scale,v2(:,1)*obj.scale*obj.framerate,'ko')
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            pbaspect([1 1 1])
            axis([0.9*min(r(:,1)) 1.1*max(r(:,1))*obj.scale 1.1*min(v2(:,1))*obj.scale*obj.framerate 1.1*max(v2(:,1))*obj.scale*obj.framerate])
            box on
            hold off
            
        end
        function plot_tracks(obj)
            %Plots the particle trajectories in time over the first image
            %of the movie.
            %Overrides plot_tracks in class particles
                        
            
            traj=get_trace(obj);
            if ~isnan(obj.im)
                picture=obj.im;
                %Assert image is 8-bit
                if length(size(picture))>2
                    %Convert to grayscale
                    picture = rgb2gray(picture);
                end
            end
              
            figure
            hold on
            box on
            if ~isnan(obj.im)
                colormap('gray'),imagesc(picture);
                axis([0 length(picture(1,:)) 0 length(picture(:,1))]);
            end
            
                trj1=traj(1:obj.NOF,1:2);
                trj2=traj(obj.NOF+1:end,1:2);
                plot(trj1(:,1),trj1(:,2),'r-','DisplayName','Particle 1: chasee','LineWidth',1);
                plot(trj2(:,1),trj2(:,2),'b-','DisplayName','Particle 2: chaser','LineWidth',1);
            
            xlabel('X (pixels)','FontSize',18)
            ylabel('Y (pixels)','FontSize',18)
            lgd.Color=[0.5 0.5 0.5];
            lgd.TextColor=[1 1 1];
            lgd=legend('show','Location','northeastoutside');
            pbaspect([1 1 1])
            %axis([800 1050 500 750])
            hold off
            
        end
        function obj=analyze_spiral(obj)
            %use only on spiral trajectories of dropletpairs
            %calculates the difference in speed between the predator and
            %prey particle
            %calculates the local curvature of the prey particle trajectory
            %plots the speed difference vs local curvature
            
            %first calculate curvature
            for t=2:obj.NOF-1
                pred1=obj;
            end
            
        end
    end %end methods
end %end class