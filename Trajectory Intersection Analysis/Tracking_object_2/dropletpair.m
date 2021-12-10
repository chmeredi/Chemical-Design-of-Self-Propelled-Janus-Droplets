classdef dropletpair<droplets
    %Droplets class is a sublcass of particles
    %Extra properties
    %   - tr_rot: rotated trajectory starting from x axis with center on origin
    %   - r_diff: list of distance between pair as function of time
    %   - r_cent: position of particle as function of time
    %   - v_rel: interparticle speed as function of interpaticle distance (interaction).
    %   - v_sum: speed of center of mass as function of interparticle distance
    %Extra methods
    %   - align_initial_frame(obj): creates tr_rot.
    %   - find_relative_displacement(obj): creates r_diff and r_cent.
    %   - find_interaction
    %   - rescale_interaction
    %   - calc_speed_vs_distance
    %Created on 04-09-17 by Pepijn Moerman
    %Last modified: -
    
    properties
        tr_rot                  % Trajectories of pair particles, rotated in such a way that 
        %in the first frame the pair lies on the x-axis with the center on the origin
        %Of type [x y t ID] with ID is 1 and 2 for particle 1 and 2.
        r_diff_rot              % Distance between particles in x, where x is along the line of their initial interaction axis
        r_cent_rot              % Average center (c1+c2)/2 of particles in x, where x is along the line of their initial interaction axis
        r_diff                  % Distance between particles
        r_cent                  % Average center (c1+c2)/2 of particles
        v_rel_rot               % Instant speed of particles with respect to each other, along x axis
        v_sum_rot               % Instant speed of center of mass, along x axis
        v_rel                   % Instant speed of particles with respect to each other
        v_sum                   % Instant speed of center of mass
        diameters               % Droplet diameters [diameter1 ID1 ; diameter2 ID2]
        rescaled_interaction
        fit_rescaled_interaction
        speedvsdist
        speedvsdistrot
    end
    methods
        function obj=dropletpair(fname,parameters,identity)
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
        function new_obj=crop_in_time_trace(obj,tbegin,tend)
            %saves only part of the trajectory between tbegin and tend
            %assumes only two particles in object
            
            if nargin<3
                error('Please provide a begin and end time.')
            end
            
            if obj.NOP~=2
                error('This function only works if you have exactly 2 particles.')
            end
            
            if tbegin<1
                tbegin=1;
            end
            if tend>obj.NOF
                tend=obj.NOF;
            end
            if tbegin>tend
                error('Please make sure the begin time is smaller than the end time.')
            end
            
            tr1=obj.tr(1:obj.NOF,:);
            tr1=tr1(tbegin:tend,:);
            tr1(:,3)=linspace(1,length(tr1(:,3)),length(tr1(:,3)))';
            tr2=obj.tr(obj.NOF+1:2*obj.NOF,:);
            tr2=tr2(tbegin:tend,:);
            tr2(:,3)=linspace(1,length(tr1(:,3)),length(tr1(:,3)))';
            
            trcombi=[tr1; tr2];
            
            new_obj=obj;
            
            new_obj.tr=trcombi;
            new_obj.NOF=length(tr1(:,1));
        end
        function obj=align_initial_frame(obj)
            %Takes the pair of droplet trajectory and aligns using the
            %inital frame in such a way that both droplets sit on the x
            %axis and the middle point is on the origin.
            if isempty(obj.tr)
                obj.tr=track_obj(obj);
                disp('Executed tracking first because trace was still empty')
            end
            tr=obj.tr;
            tr1=tr(find(tr(:,4)==1),1:4);
            tr2=tr(find(tr(:,4)==2),1:4);
            
            %Find first frame that isn't NaN
            breakval=0;
            start=1;
            for i=1:obj.NOF
                if breakval==0
                    if isnan(tr1(i,1)) || isnan(tr2(i,1))
                        start=i+1;
                    else
                        breakval=1;
                    end
                end
            end
                        
            %Find initial position of center of mass and calculate reduced
            %trajectories (substract inital position)
            pos_initial=[mean([tr1(start,1) tr2(start,1)]) mean([tr1(start,2) tr2(start,2)])];
            tr1_red=[tr1(:,1)-pos_initial(1) tr1(:,2)-pos_initial(2) tr1(:,3)];
            tr2_red=[tr2(:,1)-pos_initial(1) tr2(:,2)-pos_initial(2) tr1(:,3)];
            
            %Rotate all so that angle between x-axis and droplets is 0.
            %First find angle
            tr2_init=[tr2(start,1)-tr1(start,1) tr2(start,2)-tr1(start,2)];
            alpha=atan2(tr2_init(2),tr2_init(1));
            %Then rotate all
            tr1_rot=NaN(length(tr1(:,1))+1-start,1);
            tr2_rot=NaN(length(tr1(:,1))+1-start,1);
            for j=start:length(tr1(:,1))
                R=[cos(-alpha)  -sin(-alpha) ; sin(-alpha)  cos(-alpha)];
                tr1_rot(j,1:2)=R*tr1_red(j,1:2)';
                tr2_rot(j,1:2)=R*tr2_red(j,1:2)';
            end
            %Save in object
            tr1_rot_tot=[tr1_rot(:,1:2) tr1(:,3:4)];
            tr2_rot_tot=[tr2_rot(:,1:2) tr2(:,3:4)];
            obj.tr_rot=[tr1_rot_tot ; tr2_rot_tot];
            
            %Plot rotated trajectories
            figure
            hold on
            box on
            plot(tr1_rot(:,1),tr1_rot(:,2),'Color','r','DisplayName',['Particle ' num2str(1)]);
            plot(tr2_rot(:,1),tr2_rot(:,2),'Color','b','DisplayName',['Particle ' num2str(2)]);
            xlabel('X (pixels)','FontSize',18)
            ylabel('Y (pixels)','FontSize',18)
            %lgd.Color=[0.5 0.5 0.5];
            lgd.TextColor=[0 0 0];
            lgd=legend('show');
            hold off
        end
        function obj=find_relative_displacement_rot(obj)
            % Takes as input the aligned trajectory tr_rot.
            % Considers only motion along the direction of particle vector(x)
            % Outputs particle distance vs time and central position vs
            % time.
            
            %Calculate obj.tr_rot if it hasn't been done yet
            if isempty(obj.tr_rot)
                obj=align_initial_frame(obj);
            end
            
            %Calculate distance vs time
            tr_roti=obj.tr_rot;
            tr_rot1=tr_roti(find(tr_roti(:,4)==1),1);
            tr_rot2=tr_roti(find(tr_roti(:,4)==2),1);
            obj.r_diff_rot=[abs(tr_rot2-tr_rot1) tr_roti(find(tr_roti(:,4)==1),3)];
            
            %Calculate position vs time
            obj.r_cent_rot=[(tr_rot2+tr_rot1)/2 tr_roti(find(tr_roti(:,4)==1),3)];
            
            %Plot the values
            figure
            hold on
            box on
            
            plot(obj.r_diff_rot(:,2),obj.r_diff_rot(:,1),'ko','DisplayName','Difference');
            plot(obj.r_cent_rot(:,2),obj.r_cent_rot(:,1),'ro','DisplayName','Center');
            xlabel('t (frames)','FontSize',18)
            ylabel('X (pixels)','FontSize',18)
            lgd.Color=[0.5 0.5 0.5];
            lgd.TextColor=[1 1 1];
            lgd=legend('show');
            hold off

        end
        function obj=find_relative_displacement(obj)
            % Takes as input the aligned trajectory tr.
            % Outputs particle distance vs time and central position vs
            % time.
            
            %Calculate distance vs time
            tr_i=obj.tr;
            tr_i1x=tr_i(find(tr_i(:,4)==1),1);
            tr_i2x=tr_i(find(tr_i(:,4)==2),1);
            tr_i1y=tr_i(find(tr_i(:,4)==1),2);
            tr_i2y=tr_i(find(tr_i(:,4)==2),2);
            obj.r_diff=[sqrt((tr_i1x-tr_i2x).^2+(tr_i1y-tr_i2y).^2) tr_i(find(tr_i(:,4)==1),3)];
            
            %Calculate position vs time
            obj.r_cent=[(tr_i1x+tr_i2x)/2 (tr_i1y+tr_i2y)/2 tr_i(find(tr_i(:,4)==1),3)];
            
            %Plot the values
            figure
            hold on
            box on
            
            plot(obj.r_diff(:,2),obj.r_diff(:,1),'ko','DisplayName','Difference');
            plot(obj.r_cent(:,2),obj.r_cent(:,1),'ro','DisplayName','Center');
            xlabel('t (frames)','FontSize',18)
            ylabel('X (pixels)','FontSize',18)
            lgd.Color=[0.5 0.5 0.5];
            lgd.TextColor=[1 1 1];
            lgd=legend('show');
            hold off

        end
        function obj=find_interaction_rot(obj)
            %Takes the displacement curves tr_cent and tr_rot
            %and calculates from this the interaction potential U vs r
            %It calculates U (derivative of distance with time) vs distance
            %r as interaction potential AND U_tot (derivative of central
            %position with time) vs distance r.
            
            if isempty(obj.r_cent_rot) || isempty(obj.r_diff_rot)
                obj=find_relative_displacement_rot(obj);
            end
            
            % Obtain the reduced tracks
            r_diff_i=obj.r_diff_rot;
            r_cent_i=obj.r_cent_rot;
            
            % Calculate speeds from reduced tracks
            v_rel_i=NaN(obj.NOF,2);
            v_sum_i=NaN(obj.NOF,2);
            for i=2:obj.NOF
                v_rel_i(i,1)=r_diff_i(i,1)-r_diff_i(i-1,1);
                v_rel_i(i,2)=(r_diff_i(i,1)+r_diff_i(i-1,1))/2;
                v_sum_i(i,1)=r_cent_i(i,1)-r_cent_i(i-1,1);
                v_sum_i(i,2)=(r_diff_i(i,1)+r_diff_i(i-1,1))/2;
            end
            
            %Make sure total speed is always positive (in direction of
            %chased droplet)
            if nansum(v_sum_i(:,1))<0
                v_sum_i(:,1)=-v_sum_i(:,1);
            end
            
            %Save calculated speeds in object
            obj.v_rel_rot=v_rel_i;
            obj.v_sum_rot=v_sum_i;
            
            %Plot speeds
            figure
            hold on
            box on
            plot(v_rel_i(:,2)*obj.scale,v_rel_i(:,1)*obj.scale*obj.framerate,'ko','DisplayName','Interaction');
            plot(v_sum_i(:,2)*obj.scale,v_sum_i(:,1)*obj.scale*obj.framerate,'ro','DisplayName','Total speed');
            xlabel('r (\mum)','FontSize',18)
            ylabel('Speed (\mum/s)','FontSize',18)
            legend('show')
            %lgd.Color=[0.5 0.5 0.5];
            %lgd.TextColor=[1 1 1];
            hold off
            
            %Plot speeds vs 1/r^2
            figure
            hold on
            box on
            plot(1./v_rel_i(:,2).^2*obj.scale.^2,v_rel_i(:,1)*obj.scale*obj.framerate,'ko','DisplayName','Interaction');
            plot(1./v_sum_i(:,2).^2*obj.scale.^2,v_sum_i(:,1)*obj.scale*obj.framerate,'ro','DisplayName','Total speed');
            xlabel('1/r^2 (\mum^-^2)','FontSize',18)
            ylabel('Speed (\mum/s)','FontSize',18)
            legend('show')
            %lgd.Color=[0.5 0.5 0.5];
            %lgd.TextColor=[1 1 1];
            hold off
            
        end
        function obj=find_interaction(obj)
            %Takes the displacement curves r_diff and r_cent
            %and calculates from this the interaction potential U vs r
            %It calculates U (derivative of distance with time) vs distance
            %r as interaction potential AND U_tot (derivative of central
            %position with time) vs distance r.
            
            if isempty(obj.r_cent) || isempty(obj.r_diff)
                obj=find_relative_displacement(obj);
            end
            
            % Obtain the reduced tracks
            r_diff_i=obj.r_diff;
            r_cent_i=obj.r_cent;
            
            % Calculate speeds from reduced tracks
            v_rel_i=NaN(obj.NOF,2);
            v_sum_i=NaN(obj.NOF,2);
            for i=2:obj.NOF
                v_rel_i(i,1)=r_diff_i(i,1)-r_diff_i(i-1,1);
                v_rel_i(i,2)=(r_diff_i(i,1)+r_diff_i(i-1,1))/2;
                v_sum_i(i,1)=r_cent_i(i,1)-r_cent_i(i-1,1);
                v_sum_i(i,2)=(r_diff_i(i,1)+r_diff_i(i-1,1))/2;
            end
            
            %Make sure total speed is always positive (in direction of
            %chased droplet)
            if nansum(v_sum_i(:,1))<0
                v_sum_i(:,1)=-v_sum_i(:,1);
            end
            
            %Save calculated speeds in object
            obj.v_rel=v_rel_i;
            obj.v_sum=v_sum_i;
            
            %Plot speeds
            figure
            hold on
            box on
            plot(v_rel_i(:,2),v_rel_i(:,1),'ko','DisplayName','Interaction');
            plot(v_sum_i(:,2),v_sum_i(:,1),'ro','DisplayName','Total speed');
            xlabel('r (pixels)','FontSize',18)
            ylabel('Speed (pixels/frame)','FontSize',18)
            legend('show')
            %lgd.Color=[0.5 0.5 0.5];
            %lgd.TextColor=[1 1 1];
            hold off
            
        end
        function obj=measure_size_manually(obj)
            %Shows image to user and measures sizes of droplets
            %Saves size to object
            
            if isempty(obj.im)
                error('There is no image attached to this object, so no size can be determined.')
            end
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            obj.diameters=cell(2,2);
            
            for p=1:obj.NOP
                test='n';
                while test=='n'
                    disp('Indicate two positions on the exterior of droplet 1, opposite from each other')
                    if ishandle(1)
                        close(1)
                    end
                    figure(1)
                    hold on
                    imagesc(obj.im);
                    [x, y] = ginput(2);
                    hold off
                    diameter=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
                    center=[mean([x(2) x(1)]),mean([y(2) y(1)])];
                    close(1)
                    figure(1)
                    hold on
                    imagesc(obj.im)
                    viscircles(center,diameter/2);
                    hold off

                    message='Is the measurement good (y) or do you want to try again (n): '; 
                    acceptables=['y' 'n'];
                    test=get_textual_input(message,acceptables);
                end
                %Find which droplet is which
                distance1=sqrt((center(1)-obj.tr(1,1))^2+(center(2)-obj.tr(1,2))^2);
                distance2=sqrt((center(1)-obj.tr(obj.NOF+1,1))^2+(center(2)-obj.tr(obj.NOF+1,2))^2);
                if distance1<distance2
                    ID=1;
                else
                    ID=2;
                end
                obj.diameters{p,1}=ID;
                obj.diameters{p,2}=diameter;
            end
            close(1);
        end
        function obj=scale_attraction_data(obj)
            %Scales the attraction data by size
            
            if isempty(obj.r_cent)
                obj=find_interaction(obj);
            end
            if isempty(obj.diameters)
                obj=measure_size_manually(obj);
            end
            
            d_mean=mean([obj.diameters{1,2} obj.diameters{2,2}]);
            interaction=obj.v_rel;
            interaction(find(interaction(:,2)<1.2*d_mean),:)=NaN;
            interaction(:,1)=interaction(:,1)*d_mean;
            interaction(:,2)=interaction(:,2)/d_mean;
            
            interaction1=-interaction(:,1);
            interaction2=-interaction(:,2);
            interaction1(isnan(interaction1))=[];
            interaction2(isnan(interaction2))=[];
            
            fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',0,...
               'Upper',Inf,...
               'StartPoint',1);
            ft = fittype('a-2*x','options',fo);
            
            f1=fit(log10(interaction2),log10(interaction1),ft)
                       
            figure
            hold on
            xlabel('Log(r/d) ','FontSize',18)
            ylabel(' Log(Normalized speed) ','FontSize',18)
            plot(log10(-interaction(:,2)),log10(-interaction(:,1)),'ko');
            plot(linspace(0.8*min(log10(interaction2)),1.2*max(log10(interaction2)),100),f1.a-2*linspace(0.8*min(log10(interaction2)),1.2*max(log10(interaction2)),100),'r-');
            hold off
             
            obj.rescaled_interaction=interaction;
            obj.fit_rescaled_interaction=f1.a;
            
        end
        function obj=calc_speed_vs_distance(obj)
            %Takes the two droplets, calculates their distance and absolute
            %speed at each point and plots one against the other
            
            x1=obj.tr(1:obj.NOF,1);
            y1=obj.tr(1:obj.NOF,2);
            x2=obj.tr(obj.NOF+1:end,1);
            y2=obj.tr(obj.NOF+1:end,2);
            distance=sqrt((x2-x1).^2+(y2-y1).^2);
            
            speed1=sqrt((x1(2:end)-x1(1:end-1)).^2+(y1(2:end)-y1(1:end-1)).^2);
            speed2=sqrt((x2(2:end)-x2(1:end-1)).^2+(y2(2:end)-y2(1:end-1)).^2);
            
            figure
            hold on
            box on
            plot(distance(1:end-1),speed1,'b.');
            plot(distance(1:end-1),speed2,'r.');
            hold off
            
            obj.speedvsdist=[distance(1:end-1) speed1 speed2];
            
        end
        function obj=calc_speed_vs_distance_rot(obj)
            %Takes the two droplets, calculates their distance and absolute
            %speed at each point and plots one against the other
            
            if isempty(obj.tr_rot)
                obj=align_initial_frame(obj);
            end
            
            x1=obj.tr_rot(1:obj.NOF,1);
            y1=obj.tr_rot(1:obj.NOF,2);
            x2=obj.tr_rot(obj.NOF+1:end,1);
            y2=obj.tr_rot(obj.NOF+1:end,2);
            distance=sqrt((x2-x1).^2+(y2-y1).^2);
            
            speed1=x1(2:end)-x1(1:end-1);
            speed2=x2(2:end)-x2(1:end-1);
            
            figure
            hold on
            box on
            plot(distance(1:end-1)*obj.scale,speed1*obj.scale*obj.framerate,'bo');
            plot(distance(1:end-1)*obj.scale,speed2*obj.scale*obj.framerate,'ro');
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            hold off
            
            figure
            hold on
            box on
            loglog(distance(1:end-1),speed1,'bo');
            loglog(distance(1:end-1),speed2,'ro');
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            hold off
            
            obj.speedvsdistrot=[distance(1:end-1) speed1 speed2];
            
        end
        function obj=analyze_rot(obj)
            %combines the three rotated analysis codes for convenience
            obj=align_initial_frame(obj);
            obj=find_relative_displacement_rot(obj);
            obj=find_interaction_rot(obj);
       end
        
    end
end