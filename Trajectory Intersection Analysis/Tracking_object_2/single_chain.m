classdef single_chain < chains
    %single_chain class is a subclass of chains
    %it is used for analyzing the motion of 1D particle clusters
    %
    %Make sure this code is in the same folder as particles.m and all
    %helper codes
    %
    %Extra methods
    % - obj=single_chains()                         Initializes object
    % - obj=find_chainlength(obj)                   finds chainlength
    % - output=get_center_of_mass(obj)              returns the trajectory of the center of mass in the lab frame  
    % - obj=find_reduced_tracks(obj)                overrides function of particles. Now calculates trajectory of internal bond angles 
    % - obj=plot_reduced_tracks(obj)                plots trajectory of internal bond angle as a 2D histogram. Color codes point density. 
    % - plot_reduced_tracks_fast(obj)               plots trajectory of internal bond angle without color code 
    % - obj=find_angles(obj)                        calculates internal bond angle in each frame 
    % - plot_angle_in_time(obj)                     plot internal bond angle as function of time 
    % - plot_angle_histogram(obj)                   plots histogram of all bond angles 
    % - obj=find_center_of_mass(obj)                finds the center of mass of the cluster
    % - obj=find_angle_MSD(obj)                     calculates the angular mean squared displacement 
    % - plot_angle_MSD(obj)                         plots angular mean squared displacement 
    % - obj=find_endvector(obj)                     calculates length and direction of vector pointing from particle 1 to 3  
    % - plot_endvector_histogram(obj)               plots histogram of all end-to-end distances 
    % - obj=find_contourlength(obj)                 calculate the length of the chain from beginning to end walking through each particle 
    % - output=store_cmass(obj)                     stores the position of the center of mass in lab frame as an object of type particle 
    % - output=store_cmass_red(obj)                 stores the position of the center of mass in particle frame as type particle 
    % - plot_orient_in_time(obj)
    % - plot_orient_histogram(obj,n)                plots a histogram of all particle orientations with n bins 
    % - plot_angle_displacement_histogram(obj,n)    plots a histogram of all the angular displacements with n bins. 
    % - obj=calc_diff_tensor(obj)                   calculate diffusion tensor 
    % - obj=find_tr_mass_red(obj)                   returns the trajectory of the center of mass in the particle frame
    % - plot_scatter_correlations(obj)              makes scatterplots of all degrees of freedom versus each other 
    % - obj=crop_in_time(obj,t_min,t_max)           throws away all information between t_min and t_max 
    % - output=average_correlations(list)           returns an object containing the average diffusion tensor of all objects in list 
    % - plot_tr_mass_red(obj)
    % - displacement_histogram_red(obj)
    % - plot_tracks_plus_orient(obj)
    % - plot_image_orient(obj)
    % - plot_DxDy_trend(obj)
    % -
    %
    %Created on 19-07-17 by Pepijn Moerman
    %Last modified: 31-01-2018
    %
    properties
        chainlength     %number of particles in chain
        angles          %internal angles of chain
        orients         %external angles of the chain
        tr_reduced      %reduced trajectory of each particle in the chain with 2 neighbours
        tr_mass         %trajectory of the center of mass of the chain
        tr_mass_red     %trajectory of the center of mass, rotated to match the end-to-end vector and the x axis.
        angle_MSD       %angular mean squared displacement
        endvector       %end-to-end distance of the chain in each frame in um
        diffs           %diffusion tensor
        errors          %error bar on diffusion tensor as obtained from fit to first n points
        corrs           %correlation functions between all degrees of freedom
        contourlength   %contourlength of the chain in each frame in um
    end
    methods
        function obj=single_chain(fname,parameters,identity)
            %single_chain needs to be an object of type chains
            if nargin==1
                parameters=[];
            elseif nargin<1
                parameters=[];
                fname=[];
            end
            if nargin<3
                identity='y';
            end
            obj=obj@chains(fname,parameters,identity);
        end
        function obj=find_chainlength(obj)
            %creates chainlength indicating the number of particles in the
            %chain
            
            %if num_neighs is empty, create num_neighs first
            if isempty(obj.num_neighs)
                if isempty(obj.cutoff)
                    obj=find_neighbours(obj,40);
                    disp('The matrix num_neighs did not exist yet, so we created this matrix first. We used the default cutoff value of 40.')
                    disp('Please run find_neighbours(obj,cutoff) first with a cutoff value of your choosing to improve performance.')
                else
                    obj=find_neighbours(obj,obj.cutoff);
                    disp(['The matrix num_neighs did not exist yet, so we created this matrix first. We used the cutoff value ' num2str(obj.cutoff) ' that was stored in the object.'])
                    disp('Please run find_neighbours(obj,cutoff) first with a cutoff value of your choosing to improve performance.')
                end
            end
            
            number=obj.num_neighs;
            
            two_neighbours=find(number==2);
            one_neighbour=find(number==1);
            if ~isempty(two_neighbours)
                obj.chainlength=length(two_neighbours)+2;
            elseif ~isempty(one_neighbour)
                obj.chainlength=2;
                disp('Your chain only contains two particles');
            else
                obj.chainlength=1;
                disp('Your chain only contains one particle');
            end
        end
        function cofm=get_center_of_mass(obj)
            %ouputs the center of mass trajectory of the chain
            if isempty(obj.tr_mass)
                find_center_of_mass(obj);
            end
            cofm=obj.tr_mass;
        end
        function obj=find_reduced_tracks(obj,noshow)
            %This function gives the tracks around the frame of reference
            %of particle each central particle
            %The input variable noshow indicates whether or not to plot the
            %reduced tracks. Give (1) for yes and (0) for no. Default is
            %no.
            if nargin<2
                noshow=0;
            end
            if isempty(obj.chainlength)
                obj=find_chainlength(obj);
            end
            trace=obj.tr;
            if obj.chainlength<3
                disp('The chain is less than 3 particles long, so there is no particle with 2 neighbours to define any angle for.')
                tr_red_temp=NaN;
            else
                tr_red_temp=[];
                number=obj.num_neighs;
                two_neighbours=find(number==2);
                obj.tr_reduced={};
                for p=1:obj.NOP
                    if any(two_neighbours==p)
                        ID0=p;
                        neighbours_p=obj.neighbours{p};
                        ID1=neighbours_p(1);
                        ID2=neighbours_p(2);
                        
                        trj_0=[trace(find(trace(:,4)==ID0),1:2) trace(find(trace(:,4)==ID0),3)];
                        trj_1=[trace(find(trace(:,4)==ID1),1:2) trace(find(trace(:,4)==ID1),3)];
                        trj_2=[trace(find(trace(:,4)==ID2),1:2) trace(find(trace(:,4)==ID2),3)];
                        count=0;
                        trj_0_red=[];
                        trj_1_red=[];
                        trj_2_red=[];
                        for i=1:max([length(trj_0) length(trj_1) length(trj_2)])
                            if any(trj_0(:,3)==i) && any(trj_1(:,3)==i) && any(trj_2(:,3)==i)
                                count=count+1;
                                trj_0_red(count,1:2)=[0 0];
                                trj_1_red(count,1:2)=trj_1(find(trj_1(:,3)==i),1:2)-trj_0(find(trj_0(:,3)==i),1:2);
                                trj_2_red(count,1:2)=trj_2(find(trj_2(:,3)==i),1:2)-trj_0(find(trj_0(:,3)==i),1:2);
                            end
                        end
                        
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
                            title(['Particle ' num2str(ID0)])
                            plot(trj_0_red(:,1)*obj.scale,trj_0_red(:,2)*obj.scale,'ko')
                            plot(trj_1_rot(:,1)*obj.scale,trj_1_rot(:,2)*obj.scale,'r-')
                            plot(trj_2_rot(:,1)*obj.scale,trj_2_rot(:,2)*obj.scale,'b-')
                            hold off
                        end
                        tr_red_temp=[trj_0_red(:,1:2) ones(length(trj_0_red(:,1)),1)*ID0 ; trj_1_rot(:,1:2) ones(length(trj_1_rot(:,1)),1)*ID1 ; trj_2_rot(:,1:2) ones(length(trj_2_rot(:,1)),1)*ID2];
                        obj.tr_reduced{p}=tr_red_temp;
                    else
                        obj.tr_reduced{p}=NaN;
                    end
                    
                end
            end
        end
        function obj=plot_reduced_tracks(obj,cutoffvalue)
            %plots the reducted trajectory of each particle in the chain
            %with 2 neighbours
            %Cutoffvalue gives distance between points used to give
            %determine color scheme
            %Default cutoffvalue=2;
            if nargin==1
                cutoffvalue=2;
            end
                        
            %if tr_red exists, use object saved data, else find tr_red
            if isempty(obj.tr_reduced)
                obj=find_reduced_tracks(obj,1);
            else
                number=obj.num_neighs;
                two_neighbours=find(number==2);
                %For each particle
                figure
                count1=floor(sqrt(obj.chainlength-2))+1;
                count2=floor(sqrt(obj.chainlength-2))+1;
                count=0;
                for p=1:obj.NOP
                    %But only if the particle has 2 neighbours
                    if any(two_neighbours==p)
                        count=count+1;
                        %Obtain the reduced trajectories
                        tr_red_p=obj.tr_reduced{p};
                        neighbours_p=obj.neighbours{p};
                        ID0=p;
                        ID1=neighbours_p(1);
                        ID2=neighbours_p(2);
                        trj_0_red=tr_red_p(find(tr_red_p(:,3)==ID0),1:2);
                        trj_1_red=tr_red_p(find(tr_red_p(:,3)==ID1),1:2);
                        trj_2_red=tr_red_p(find(tr_red_p(:,3)==ID2),1:3);
                        
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
                        cmap=jet(max_num);
                        
                        %Plot the reduced trajectories
                        subplot(count1,count2,count)
                        hold on
                        box on
                        edge=obj.cutoff*obj.scale;
                        axis([-edge edge -edge edge])
                        xlabel('X (\mum)')
                        ylabel('Y (\mum)')
                        title(['Particle ' num2str(ID0)])
                        axis([-1.1*max(abs(trj_2_red(:,1)))*obj.scale 1.1*max(abs(trj_2_red(:,1)))*obj.scale -1.1*max(abs(trj_2_red(:,2)))*obj.scale 1.1*max(abs(trj_2_red(:,2)))*obj.scale]);
                        %plot(trj_0_red(:,1)*obj.scale,trj_0_red(:,2)*obj.scale,'ko')
                        %plot(trj_1_red(:,1)*obj.scale,trj_1_red(:,2)*obj.scale,'r-')
                        %plot the last reduced trajectory in colors
                        %depending on the number of neighbours
                        line([0 0], [-100 100],'Color','k');
                        line([-100 100], [0 0],'Color','k');
                        line([0 -100],[0 -173.2],'Color','r');
                        line([0 -100],[0 173.2],'Color','r');
                        for i=1:length(trj_2_red(:,1))
                            num_neigh_i=round(num_of_neighs(i));
                            if num_neigh_i<1
                                num_neigh_i=1;
                            end
                            plot(trj_2_red(i,1)*obj.scale,trj_2_red(i,2)*obj.scale,'o','Color',cmap(num_neigh_i,:))
                            %plot(trj_2_red(i,1)*obj.scale,trj_2_red(i,2)*obj.scale,'o')
                        end
                        axis([-edge edge -edge edge])
                        hold off
                    end
                end
            end
        end
        function plot_reduced_tracks_fast(obj)
            %plots the reducted trajectory of each particle in the chain
            %with 2 neighbours
            %no colorscheme to expedite calculation
            
            if isempty(obj.tr_reduced)
                obj=find_reduced_tracks(obj);
            end
            
            %find particles with two neighbours
            number=obj.num_neighs;
            two_neighbours=find(number==2);
            %for all particles with two neighbours
            for p=1:obj.NOP
                if any(two_neighbours==p)
                    %Obtain the reduced trajectories
                    tr_red_p=obj.tr_reduced{p};
                    neighbours_p=obj.neighbours{p};
                    ID0=p;
                    ID2=neighbours_p(2);
                    trj_2_red=tr_red_p(find(tr_red_p(:,3)==ID2),1:3);
                    
                    %Plot the reduced trajectories
                    figure
                    hold on
                    box on
                    axis([-1.1*max(trj_2_red(:,1))*obj.scale 1.1*max(trj_2_red(:,1))*obj.scale -1.1*max(abs(trj_2_red(:,2)))*obj.scale 1.1*max(abs(trj_2_red(:,2)))*obj.scale])
                    xlabel('X (\mum)')
                    ylabel('Y (\mum)')
                    title(['Particle ' num2str(ID0)])
                    plot(trj_2_red(:,1)*obj.scale,trj_2_red(:,2)*obj.scale,'b-')
                    %plot cross around zero
                    line([0 0], [-100 100],'Color','k');
                    line([-100 100], [0 0],'Color','k');
                    hold off
                end
            end

        end
        function obj=find_angles(obj)
            %calculates the internal angle of each colloid with two
            %neighbours
            if isempty(obj.tr_reduced)
                obj=find_reduced_tracks(obj);
            end
            if obj.chainlength<3
                disp('The chain is less than 3 particles long, so there is no particle with 2 neighbours to define any angle for.')
                obj.angles=NaN;
            else
                obj.angles=cell(obj.chainlength-2,2);
                count=0;
                for p=1:obj.NOP
                    if obj.num_neighs(p)==2
                        ID=p;
                        tr_red_temp=obj.tr_reduced{1,p};
                        tr_red_p=tr_red_temp(2*obj.NOF+1:3*obj.NOF,1:2);
                        angles_p=atan2(tr_red_p(:,2),tr_red_p(:,1));
                        count=count+1;
                        obj.angles{count,1}=ID;
                        obj.angles{count,2}=angles_p;
                    end
                end
            end
        end
        function plot_angle_in_time(obj)
            %Plots a graph for each particle with two neighbours
            %Graph contains the angle in time
            if ~iscell(obj.angles)
                obj=find_angles(obj);
            end
            for p=1:obj.chainlength-2
                anglep=obj.angles{p,2};
                figure
                hold on
                title(['Particle ' num2str(p)])
                plot(linspace(1,length(anglep),length(anglep))/obj.framerate,anglep/pi,'k-');
                xlabel('Time (s)')
                ylabel('Internal angle (\pi radians)')
                yticks([-1 -0.5 0 0.5 1]);
                box on
                hold off
            end
            
        end
        function plot_angle_histogram(obj,bins,individuals)
            %Plots a histogram for each particle with two neighbours
            %Histogram of all accessed angles
            %individuals indicates whether (1) or not (0) to show the
            %histogram of individual bond angles in the case of more than 1
            %bond. Default is no show (0)
            
            if nargin<3
                individuals=0;
            end
            
            if ~iscell(obj.angles)
                obj=find_angles(obj);
            end
            
            angle_tot=[];
            figure(10)
            count1=floor(sqrt(obj.chainlength-2))+1;
            count2=floor(sqrt(obj.chainlength-2))+1;
            for p=1:obj.chainlength-2
                anglep=obj.angles{p,2};
                if individuals~=0
                    subplot(count1,count2,p)
                    hold on
                    title(['Particle ' num2str(p)])
                    xlabel('angle (\pi radians)')
                    ylabel('P(\alpha)')
                    xticks([-1/2 0 1/2]);
                    xticklabels({'-1/2','0','1/2'});
                    axis([-2/3 2/3 -1 1])
                    axis 'auto y'
                    if nargin==1
                        histogram(anglep/pi);
                    else
                        histogram(anglep/pi,bins);
                    end
                    box on
                    hold off
                end
                angle_tot=[angle_tot ; anglep];
            end
            
            if individuals==0
                close 10
            end
            
            figure
            hold on
            title('Average')
            xlabel('Internal angle (\pi radians)')
            if nargin==1
                histogram(angle_tot/pi,'FaceColor','b','FaceAlpha',1);
            else
                histogram(angle_tot/pi,bins,'FaceColor','b','FaceAlpha',1);
            end
            xticks([-1 -2/3 -1/2 -1/3 -1/6 0 1/6 1/3 1/2 2/3 1]);
            xticklabels({'-1', '-2/3','-1/2','-1/3','-1/6','0','1/6','1/3','1/2','2/3', '1'});
            axis([-1 1 -1 1]);
            axis 'auto y'
            ylabel('P(\alpha)')
            box on
            hold off
            
        end
        function obj=find_center_of_mass(obj)
            %Takes the particle chain coordinates and finds the center of
            %mass at each frame. Returns the trajectory of the center of
            %mass
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            trace=get_trace(obj);
            nframes=obj.NOF;
            trace_mass=NaN(nframes,4);
            for i=1:nframes
                trj=trace(find(trace(:,3)==i),1:2);
                x_mean=mean(trj(:,1));
                y_mean=mean(trj(:,2));
                trace_mass(i,1:4)=[x_mean y_mean i 1];
            end
            obj.tr_mass=trace_mass;
        end
        function obj=find_angle_MSD(obj)
            %Considers the interaction angle between three particles as a
            %random walk and calculates the MSD of the motion
            if isempty(obj.angles)
                obj=find_angles(obj);
            end
            nparts=obj.chainlength-2;
            nframes=obj.NOF;
            framerate_temp=obj.framerate;
            obj.angle_MSD=cell(nparts+1,2);
            MSD_tot=NaN(floor(nframes/3)+1,nparts);
            for i=1:nparts
                % Select the part of the file about particle i
                trj = obj.angles{i,2};
                if ~isempty(trj)
                    
                    % Find thelength of the trajectory of particle i
                    tmax=length(trj(:,1));
                    
                    % Define empty matrix for MSD
                    MSD_temp=NaN(floor(tmax/3)+1,nparts);
                    
                    % For a time interval of 1 to 1/3 of the length of the trajectory
                    for t=0:floor(tmax/3)
                        MSD_all=NaN(floor(tmax/3)+1,2);
                        %Create two matrixes with a difference t between them
                        a_frontstep=trj(1+t:length(trj),1);
                        a_backstep=trj(1:length(trj)-t,1);
                        %Substract to obtain MSD
                        MSD_all=(a_frontstep(:,1)-a_backstep(:,1)).^2;
                        %Average to obtain average MSD for t and save
                        MSD_temp(t+1,1)=t/framerate_temp;
                        MSD_temp(t+1,2)=nanmean(MSD_all(:,1));
                    end
                    
                end
                %End loop over particle number
                obj.angle_MSD{1+i,1}=i;
                obj.angle_MSD{1+i,2}=MSD_temp(:,1:2);
                size(MSD_temp)
                size(MSD_tot)
                MSD_tot(:,i)=MSD_temp(:,2);
            end
            obj.angle_MSD{1,1}='Average';
            obj.angle_MSD{1,2}=[MSD_temp(:,1) mean(MSD_tot(:,:),2)];
        end
        function plot_angle_MSD(obj)
            %Plots the angular MSD on lin lin and log log plot. If more
            %than one particle has defined angles, it plots the average and
            %the individual MSD's separately
            if isempty(obj.angle_MSD)
                obj=find_angle_MSD(obj);
            end
            if obj.chainlength>3
                %Plot the individual MSDs in one figure
                figure
                cmap=jet(obj.NOP);
                hold on
                box on
                lgd=legend('show');
                for i=1:obj.chainlength-2
                    MSDi=obj.angle_MSD{i+1,2};
                    particleID=obj.angle_MSD{i+1,1};
                    plot(MSDi(:,1),MSDi(:,2),'o','Color',cmap(i,:),'DisplayName',['Particle ' num2str(particleID)])
                end
                xlabel('\Deltat (s)','FontSize',18)
                ylabel('MSAD (\mum^2/s)','FontSize',18)
                lgd.Color=[0.5 0.5 0.5];
                lgd.TextColor=[1 1 1];
                hold off
            end
            
            %Plot the average MSD in one figure
            figure
            MSD_average=obj.angle_MSD{1,2};
            hold on
            box on
            lgd=legend('show');
            plot(MSD_average(:,1),MSD_average(:,2),'ko','DisplayName','Average')
            xlabel('\Deltat (s)','FontSize',18)
            ylabel('MSD (\mum^2/s)','FontSize',18)
            lgd.Color=[0.5 0.5 0.5];
            lgd.TextColor=[1 1 1];
            hold off
            
            %Plot the aveage MSD in log log representation
            figure
            box on
            loglog(MSD_average(:,1),MSD_average(:,2),'ko')
            xlabel('\Deltat (s)','FontSize',18)
            ylabel('MSD (\mum^2/s)','FontSize',18)
        end
        function obj=find_endvector(obj)
            %Calculates the end-to-end-distance of the chain of particles
            %in each frame. obj.endvector is organized as [r, theta, t]
            if isempty(obj.chainlength)
                obj=find_chainlength(obj);
            end
            %Find the ID's of the extreme particles
            extreme=NaN(2,1);
            count=0;
            branched=0;
            for p=1:obj.chainlength
                if length(obj.neighbours{p})==1 && count<2
                    count=count+1;
                    extreme(count)=p;
                elseif obj.neighbours{p}==1
                    branched=1;
                end
            end
            
            if any(isnan(extreme))
                disp('The chain has less than 2 particles with only 1 neighbour. This means the chain is circular or the choice of cutoff was too high.')
                disp('Fix the cutoff value before calculating the end to end distance.')
                obj.endvector=NaN;
            elseif branched==1
                disp('The chain has more than 2 particles with only 1 neighbour. This means it is branched so the end-to-end-vector has no meaning.')
                obj.endvector=NaN;
            else
                trace=obj.tr;
                endvector_temp=NaN(obj.NOF,3);
                trj1=trace(find(trace(:,4)==extreme(1)),1:3);
                trj2=trace(find(trace(:,4)==extreme(2)),1:3);
                endvector_temp(:,1)=sqrt((trj2(:,1)-trj1(:,1)).^2+(trj2(:,2)-trj1(:,2)).^2);
                endvector_temp(:,2)=atan2(trj2(:,2)-trj1(:,2),trj2(:,1)-trj1(:,1));
                %endvector_temp(find(endvector_temp(:,2)<0),2)=endvector_temp(find(endvector_temp(:,2)<0),2)+pi;
                endvector_temp(:,3)=trj2(:,3);
                obj.endvector=endvector_temp;
            end
        end
        function obj=get_orients(obj,noshow)
            %calculate the unbounded angle of the chain with the labframe
            %unbounded means angle can increase over 2pi if more than one
            %rotation is done
            %See Edmond et al. PNAS 209, 17891-17896 (2012) for details.
            %produces the property orients
            
            if nargin==1
                noshow=0;
            end
            
            %check if endvector is calculated
            if isempty(obj.endvector)
                obj=find_endvector(obj);
            end
            
            %calculate the difference in angle between each timestep
            theta1=obj.endvector(2:end,2);
            theta2=obj.endvector(1:end-1,2);
            theta3=theta1-theta2;
            %if the particle makes a complete revolution, make sure theta3
            %does not overshoot to 2pi or -2pi
            turns=find(theta3(:,1)>1.8*pi);
            theta3(turns,1)=theta3(turns,1)-2*pi;
            turns=find(theta3(:,1)<-1.8*pi);
            theta3(turns,1)=theta3(turns,1)+2*pi;
            
            orient_temp=NaN(length(obj.endvector(:,1)),2);
            orient_temp(1,1)=1;
            orient_temp(1,2)=0;
            for t=1:length(obj.endvector(:,1))-1
                orient_temp(t+1,1)=t+1;
                orient_temp(t+1,2)=sum(theta3(1:t,1));
            end
            obj.orients=orient_temp;
            if noshow==0
                plot_orient_in_time(obj);
            end
            
        end
        function obj=plot_orient_in_time(obj)
            %plots the unbounded orientation of the particles in time
            
            if isempty(obj.orients)
                obj=get_orients(obj,1);
            end
            
            figure
            hold on
            box on
            axis([0 max(obj.orients(:,1))/obj.framerate -1.1*max(abs(obj.orients(:,2))) 1.1*max(abs(obj.orients(:,2)))])
            plot(obj.orients(:,1)/obj.framerate, obj.orients(:,2),'k.');
            xlabel('t (s)');
            ylabel('\theta (rad)');
            hold off
            
        end
        function plot_endvector_histogram(obj)
            %Plots a histogram of the end-to-enddistance of the chain
            if isempty(obj.endvector)
                obj=find_endvector(obj);
            end
            if isempty(obj.contourlength)
                obj=find_contourlength(obj);
            end
            figure
            hold on
            box on
            legend('show')
            histogram(obj.endvector(:,1)*obj.scale,'DisplayName','End-to-end distance')
            histogram(obj.contourlength(:,1)*obj.scale,'DisplayName','Contour length')
            xlabel('length (\mum)','FontSize',18)
            ylabel('Number','FontSize',18)
            hold off
        end
        function obj=find_contourlength(obj)
            %Finds the contour lengths of a single chain in each frame
            if isempty(obj.neighbours)
                if isempty(obj.cutoff)
                    obj=find_neighbours(obj,40);
                    disp('obj.neighbours was still empty, so we found neighbours using default cutoff of 40. First find neighbours with a cutoff of your choosing to improve performance');
                else
                    obj=find_neighbours(obj,obj.cutoff);
                end
            end
            endpoints=find(obj.num_neighs==1);
            IDbegin=endpoints(1);
            sequence=findchain_class(obj.neighbours,obj.num_neighs,IDbegin);
            dist=NaN(obj.NOF,obj.NOP-1);
            for i=1:obj.NOP-1
                IDi=sequence(i);
                IDj=sequence(i+1);
                tri=obj.tr(find(obj.tr(:,4)==IDi),1:2);
                trj=obj.tr(find(obj.tr(:,4)==IDj),1:2);
                dist(:,i)=sqrt((trj(:,1)-tri(:,1)).^2+(trj(:,2)-tri(:,2)).^2);
            end
            obj.contourlength=sum(dist,2);
        end
        function output=store_cmass(obj)
            %takes the center of mass related to this object and stores it
            %as an object of type 'particles'
            
            if isempty(obj.tr_mass)
                obj=find_center_of_mass(obj);
            end
            cmass=get_center_of_mass(obj);
            parameters=get_parameters(obj);
            fname=[obj.foldername obj.filename '.tif'];
            output_temp=particles(fname,parameters);
            output_temp=add_trace(output_temp,cmass);
            output=output_temp;
        end
        function output=store_cmass_red(obj)
            %takes the center of mass related to this object and stores it
            %as an object of type 'particles'
            
            if isempty(obj.tr_mass_red)
                obj=find_tr_mass_red(obj);
            end
            cmass=obj.tr_mass_red;
            cmass(:,4)=1;
            parameters=get_parameters(obj);
            fname=[obj.foldername obj.filename '.tif'];
            output_temp=particles(fname,parameters);
            output_temp=add_trace(output_temp,cmass);
            output=output_temp;
        end
        function plot_orient_histogram(obj,n)
            %plots a histogram of the change in orientation of the external
            %angle with the laboratory frame
            %n is the number of bins in the histogram
            
            if nargin==1
                n=100;
            end
            
            %check if endvector is calculated
            if isempty(obj.endvector)
                obj=find_endvector(obj);
            end
            
            %calculate the difference in angle between each timestep
            theta1=obj.endvector(2:end,2);
            theta2=obj.endvector(1:end-1,2);
            theta3=theta1-theta2;
            %if the particle makes a complete revolution, make sure theta3
            %does not overshoot to 2pi or -2pi
            turns=find(theta3(:,1)>1.8*pi);
            theta3(turns,1)=theta3(turns,1)-2*pi;
            turns=find(theta3(:,1)<-1.8*pi);
            theta3(turns,1)=theta3(turns,1)+2*pi;
            
            figure
            hold on
            box on
            %histogram(theta3, 'FaceColor', 'b','FaceAlpha',1)
            h=histfit(theta3(:,1),n);
            h(1).FaceColor='b';
            h(1).FaceAlpha=1;
            h(2).Color='r';
            xlabel('\Delta\theta (rad)');
            ylabel('P(\Delta\theta)');
            hold off
            
            figure
            hold on
            box on
            %histogram(theta3, 'FaceColor', 'b','FaceAlpha',1)
            h=histfit(theta3(:,1),n);
            h(1).FaceColor='b';
            h(1).FaceAlpha=1;
            h(2).Color='r';
            set(gca,'YScale','log')
            xlabel('\Delta\theta (rad)');
            ylabel('P(\Delta\theta)');
            hold off
        end
        function plot_orient_MSD(obj)
            %plot the MSD of the angle of the cluster with the external
            %lab frame
            
            if isempty(obj.corrs)
                obj=calc_diff_tensor(obj);
            end
            
            Crot=obj.corrs{3,3};
            Drot=obj.diffs{3,3};
            figure
            hold on
            box on
            ylabel('C(t) (rad^2)');
            xlabel('dt (s)');
            plot(Crot(:,1)/obj.framerate,Crot(:,2),'ko');
            plot(Crot(:,1)/obj.framerate,Crot(:,1)/obj.framerate*Drot,'r-');
            hold off
            
            
        end
        function plot_orient_histogram_continuous(obj,n)
            %plots a histogram of the change in orientation of the external
            %angle with the laboratory frame
            %uses the continuous angular displacement orient
            %n is the number of bins in the histogram
            
            if nargin==1
                n=100;
            end
            
            %check if endvector is calculated
            if isempty(obj.orients)
                obj=get_orients(obj);
            end
            
            %calculate the difference in angle between each timestep
            theta1=obj.orients(2:end,2);
            theta2=obj.orients(1:end-1,2);
            theta3=theta1-theta2;
            %if the particle makes a complete revolution, make sure theta3
            %does not overshoot to 2pi or -2pi
            turns=find(theta3(:,1)>1.8*pi);
            theta3(turns,1)=theta3(turns,1)-2*pi;
            turns=find(theta3(:,1)<-1.8*pi);
            theta3(turns,1)=theta3(turns,1)+2*pi;
            
            figure
            hold on
            box on
            %histogram(theta3, 'FaceColor', 'b','FaceAlpha',1)
            h=histfit(theta3(:,1),n);
            h(1).FaceColor='b';
            h(1).FaceAlpha=1;
            h(2).Color='r';
            xlabel('\Delta\theta (rad)');
            ylabel('P(\Delta\theta)');
            hold off
            
            figure
            hold on
            box on
            %histogram(theta3, 'FaceColor', 'b','FaceAlpha',1)
            h=histfit(theta3(:,1),n);
            h(1).FaceColor='b';
            h(1).FaceAlpha=1;
            h(2).Color='r';
            set(gca,'YScale','log')
            xlabel('\Delta\theta (rad)');
            ylabel('P(\Delta\theta)');
            hold off
        end
        function plot_angle_displacement_histogram(obj,n)
            %plots a histogram of the angular displacement of the object
            %n is the number of bins in the histogram
            
            if nargin==1
                n=100;
            end
            
            if isempty(obj.angles)
                obj=find_angles(obj);
            end
            
            
            for p=1:length(obj.angles(:,1))
                angle_temp=obj.angles{p,2};
                angle1=angle_temp(2:end,1);
                angle2=angle_temp(1:end-1,1);
                angle3=angle1-angle2;
                figure
                hold on
                box on
                %histogram(angle3(:,1),'FaceColor','b');
                h=histfit(angle3(:,1),n);
                h(1).FaceColor='b';
                h(1).FaceAlpha=1;
                h(2).Color='r';
                xlabel('\Delta\alpha (rad)');
                ylabel('P(\Delta\alpha)');
                hold off
                
                figure
                hold on
                box on
                %histogram(angle3(:,1),'FaceColor','b','FaceAlpha',1);
                h=histfit(angle3(:,1),n);
                h(1).FaceColor='b';
                h(1).FaceAlpha=1;
                h(2).Color='r';
                set(gca,'YScale','log')
                xlabel('\Delta\alpha (rad)');
                ylabel('P(\Delta\alpha)');
                hold off
                
           end
            
            
        end
        function obj=calc_diff_tensor(obj,nfit,imagesupress)
            %Cross-correlates the all degrees of freedom (chainlength + 1) 
            %of a floppy chain in 2D with each other: 
            %tx, ty, r, alpha1, alpha2... 
            %Here tx and ty are translational motion of the instantaneous 
            %center of hydrodynamic friction, r is the rotational diffusion
            %of the particle frame with respect to the lab frame and
            %alpha_i are the internal angles and the number of alphas is 
            %chainlength-2.
            %Saves slopes at small time intervals in the diffusion tensor
            %Number of points in the correlation function used to calculate
            %slope is nfit. The default value is 10.
            %Imagesupress indicates whether or not to show the correlation
            %graphs. (1) is supress, (0) is don't supress. Default is (0);
            
            %use default values when imput is not provided
            if nargin<2 || nfit<2
                nfit=2;
            end
            if nargin<3
                imagesupress=0;
            end
            
            %Check if all required data are calculated
            if isempty(obj.tr_mass_red)
                obj=find_tr_mass_red(obj); %get center of hydrodynamic friction at each time point
            end
            if isempty(obj.angles)
                obj=find_angles(obj); %get internal bond angles at each time point
            end
            if isempty(obj.orients)
                obj=get_orients(obj,1); %get orientation of body frame with respect to lab frame, at each time point
            end
            
            %Create empty cell with length = degrees of freedom (DOF)
            funcs=cell(obj.chainlength+1,1);
            %Save x position, y position and orientation in funcs
            funcs{1}=[obj.tr_mass_red(:,3) obj.tr_mass_red(:,1)];
            funcs{2}=[obj.tr_mass_red(:,3) obj.tr_mass_red(:,2)];
            funcs{3}=[obj.orients(:,1) obj.orients(:,2)];
            %Save internal angle at each time in funcs
            %Number of internal angles depends on chain length
            for i=1:obj.chainlength-2
                funcs{i+3}=[obj.endvector(:,3) obj.angles{i,2}+pi;];
            end
            
            %create an empty correlation matrix corr with size DOFxDOF
            %(DOF=degree of freedom = chainlength+1)
            corr=cell(obj.chainlength+1,obj.chainlength+1);
            %create an empty diffusion tensor of size DOFxDOF
            Diffs=cell(obj.chainlength+1,obj.chainlength+1);
            %create an empty error in diffusion tensor of size DOFxDOF
            Errors=cell(obj.chainlength+1,obj.chainlength+1);
            %for all DOF, fill correlation matrix and diff tensor
            %and plot the correlation functions
            
            %Figure preamble
            f = figure(10);
            p = uipanel('Parent',f,'BorderType','none','BackgroundColor','white');
            p.Title = 'Diffusion tensor';
            p.TitlePosition = 'centertop';
            p.FontSize = 12;
            p.FontWeight = 'bold';
            
            %Loop over all DOFs
            for i=1:obj.chainlength+1
                %Loop over all DOFs again
                for j=1:obj.chainlength+1
                    %Cross (or auto) correlate the DOFs
                    Cij=correlate_functions(funcs{i},funcs{j});
                    %Save correlation function in corr
                    corr{i,j}=Cij;
                    %Calculate slope of first nfit points
                    Dijall=(Cij(1:nfit,2)-Cij(1,2))./(Cij(1:nfit,1)-Cij(1,1));
                    %Average the slope
                    Dij=nanmean(Dijall(:,1));
                    %Get standard deviation in the slope
                    Dijerror=nanstd(Dijall(:,1));
                    %Multiply value with scale^2 and framerate if in the
                    %top left quadrant of tensor
                    if (i==1 || i==2) && (j==1 || j==2)
                        Diffs{i,j}=Dij*obj.scale^2*obj.framerate;
                        Errors{i,j}=Dijerror*obj.scale^2*obj.framerate;
                    %Multiply with scale and framerate if in the top right
                    %or bottom left quadrant of the tensor
                    elseif i==1 || i==2 || j==1 || j==2
                        Diffs{i,j}=Dij*obj.scale*obj.framerate;
                        Errors{i,j}=Dijerror*obj.scale*obj.framerate;
                    %Multiply with framerate if in the bottom right
                    %quadrant of the diffusion tensor
                    else
                        Diffs{i,j}=Dij*obj.framerate;
                        Errors{i,j}=Dijerror*obj.framerate;
                    end
                    %Plot the correlation functions
                    if imagesupress==0
                        hold on
                        box on
                        axi=subplot(obj.chainlength+1,obj.chainlength+1,(j-1)*(obj.chainlength+1)+i,'Parent',p);
                        title(['C' num2str(i) num2str(j)])
                        xlabel('dt (s)')
                        if (i==1 || i==2) && (j==1 || j==2)
                            ylabel('C(t) (\mum^2)');
                            plot(axi,Cij(:,1)/obj.framerate,Cij(:,2)*obj.scale^2,'k.');
                        elseif i==1 || i==2 || j==1 || j==2
                            ylabel('C(t) (\mum.rad)');
                            plot(axi,Cij(:,1)/obj.framerate,Cij(:,2)*obj.scale,'k.');
                        else
                            ylabel('C(t) (rad^2)');
                            plot(axi,Cij(:,1)/obj.framerate,Cij(:,2),'k.');
                        end
                        hold off
                    end
                end
            end
            if imagesupress~=0
                close 10
            end
            %save diffusion tensor and correllation function in object
            obj.diffs=Diffs;
            obj.corrs=corr;
            obj.errors=Errors;
        end
        function plot_dxdy_MSD(obj)
            %plots the diffusion in x and in y as function of time
            if isempty(obj.corrs)
                obj=calc_diff_tensor(obj);
            end
            
            Dx=obj.corrs{1,1};
            Dy=obj.corrs{2,2};
            Dxx=obj.diffs{1,1};
            Dyy=obj.diffs{2,2};
            figure
            hold on
            plot(Dx(:,1)/obj.framerate,Dx(:,2)*obj.scale^2,'bo');
            plot(Dx(:,1)/obj.framerate,Dx(:,1)/obj.framerate*Dxx,'b-')
            plot(Dy(:,1)/obj.framerate,Dy(:,2)*obj.scale^2,'ro');
            plot(Dy(:,1)/obj.framerate,Dy(:,1)/obj.framerate*Dyy,'r-')
            box on
            hold off
        end
        function obj=find_tr_mass_red(obj)
            %Calculates the reduced trajectory of the center of mass;
            %This trajectory has the end-to-end vector (long axis of the chain) point in the
            %direction of the x-axis
            if isempty(obj.endvector)
                obj=find_endvector(obj);
            end
            if isempty(obj.tr_mass)
                obj=find_center_of_mass(obj);
            end
            
            nframes=obj.NOF;
            tr_temp=NaN(nframes,3);
            tr_temp(1,1:3)=[0 0 1];
            for i=2:obj.NOF
                angle=-(obj.endvector(i-1,2)+obj.endvector(i,2))/2;
                R=[cos(angle) -sin(angle) ; sin(angle) cos(angle)];
                dx=R*(obj.tr_mass(i,1:2)-obj.tr_mass(i-1,1:2))';
                tr_temp(i,1:2)=tr_temp(i-1,1:2)+dx';
                tr_temp(i,3)=i;
            end
            obj.tr_mass_red=tr_temp;
            
            figure
            hold on
            box on
            xlabel('x (\mum)')
            ylabel('y (\mum)')
            plot(tr_temp(:,1),tr_temp(:,2),'b-')
            axis equal
%             for t=1:10:length(tr_temp(:,1))
%                 a=obj.NOP*obj.scale;
%                 b=1*obj.scale;
%                 x0=tr_temp(t,1);
%                 y0=tr_temp(t,2);
%                 w=-pi:0.01:pi;
%                 x=x0+a*cos(w);
%                 y=y0+b*sin(w);
%                 plot(x,y,'r-','LineWidth',2)
%             end
            hold off
            
        end
        function plot_tr_mass_red(obj,dt,tmin,tmax)
            %plots the reduced trajectory of the center of mass
            %reduced means x and y are in the body frame not in the lab
            %frame where x is the long axis and y the short axis.
            
            if isempty(obj.tr_mass_red)
                obj=find_tr_mass_red(obj);
            end
            
            if nargin<2
                dt=10;
            end
            if nargin<3 || tmin<2
                tmin=2;
            end
            if nargin<4 || tmax>length(obj.tr_mass_red(:,1))
                tmax=length(obj.tr_mass_red(:,1));
                disp(['Using tmax = ' num2str(tmax) '.']);
            end
                        
            
            figure
            hold on
            box on
            xlabel('x (\mum)')
            ylabel('y (\mum)')
            plot(obj.tr_mass_red(tmin:tmax,1),obj.tr_mass_red(tmin:tmax,2),'b-')
            axis equal
            L=obj.NOP*0.5*obj.scale;
            for t=tmin:dt:tmax
                %calculate instantaneous orientation of particle
                theta_t=obj.endvector(t,2);
                x0=obj.tr_mass_red(t,1);
                y0=obj.tr_mass_red(t,2);
                x1=x0+L*cos(theta_t);
                x2=x0-L*cos(theta_t);
                y1=y0+L*sin(theta_t);
                y2=y0-L*sin(theta_t);
                x=linspace(x1,x2,10);
                %calculate alpha and rotate orientation
                angle=-(obj.endvector(t-1,2)+obj.endvector(t,2))/2;
                R=[cos(angle) -sin(angle) ; sin(angle) cos(angle)];
                r1n = R*[x1-x0;y1-y0]+[x0;y0];
                r2n = R*[x2-x0;y2-y0]+[x0;y0];
                
                y=r2n(2)+(r1n(2)-r2n(2))/(r1n(1)-r2n(1))*(x-r2n(1));                
                plot(x,y,'r-','LineWidth',1)
            end
            xmin=min(obj.tr_mass_red(tmin:tmax,1));
            xmax=max(obj.tr_mass_red(tmin:tmax,1));
            ymin=min(obj.tr_mass_red(tmin:tmax,2));
            ymax=max(obj.tr_mass_red(tmin:tmax,2));
            axis([xmin xmax ymin ymax])
            pbaspect([ 1 1 1]);
            hold off
            
            
        end
        function displacement_histogram_red(obj,nbins,show_individuals)
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
            if isempty(obj.tr_mass_red)
                obj=find_tr_mass_red(obj);
            end
            trace=obj.tr_mass_red;
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
            
            trj=trace(:,1:2);
            x1=trj(1:length(trj(:,1))-1,1);
            x2=trj(2:length(trj(:,1)),1);
            dx=x2-x1;
            dx_tot=[dx_tot ; dx];
            if show_individuals==1
                histogram(dx,'FaceColor',cmap(p,:),'DisplayName',['Particle ' num2str(p)])
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
            
            trj=trace(:,1:2);
            y1=trj(1:length(trj(:,1))-1,2);
            y2=trj(2:length(trj(:,1)),2);
            dy=y2-y1;
            dy_tot=[dy_tot ; dy];
            if show_individuals==1
                histogram(dy*obj.scale,'FaceColor',cmap(p,:),'DisplayName',['Particle ' num2str(p)])
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
            h=histfit(dx_tot*obj.scale,nbins+12);
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
            h=histfit(dx_tot*obj.scale,nbins+12);
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
        function plot_tracks_plus_orient(obj,dt,tmin,tmax)
            %plots the trajectory of an object with its instantaneous
            %orientation every dt timestep away
            %default is dt=10;
            
            if isempty(obj.endvector)
                obj=find_endvector(obj);
            end
            if isempty(obj.tr_mass)
                obj=find_center_of_mass(obj);
            end
            
            if nargin<2
                dt=10;
            end
            if nargin<3 || tmin<1
                tmin=1;
            end
            if nargin<4 || tmax>length(obj.tr_mass(:,1))
                tmax=length(obj.tr_mass(:,1));
                disp(['Using tmax = ' num2str(tmax) '.']);
            end
                        
            
            figure
            hold on
            box on
            xlabel('x (\mum)')
            ylabel('y (\mum)')
            plot(obj.tr_mass(tmin:tmax,1),obj.tr_mass(tmin:tmax,2),'b-')
            axis equal
            L=obj.NOP*10*obj.scale;
            for t=tmin:dt:tmax
                theta_t=obj.endvector(t,2);
                x0=obj.tr_mass(t,1);
                y0=obj.tr_mass(t,2);
                x1=x0+L*cos(theta_t);
                x2=x0-L*cos(theta_t);
                y1=y0+L*sin(theta_t);
                y2=y0-L*sin(theta_t);
                x=linspace(x1,x2,10);
                y=y2+(y1-y2)/(x1-x2)*(x-x2);                
                plot(x,y,'r-','LineWidth',1)
            end
            hold off
        end
        function plot_image_orient(obj,frame,vectorsize,linewidth)
            %Plots framenumber 'frame' of the associated tiffstack and
            %overlays the found particles in this image
            %vectorsize is the size of the particle in micron
            %Default frame is 1;
            %Default vectorsize is 3;
            %Default linewidth is 2;
            
            if nargin==1
                frame=1;
                vectorsize=2;
                linewidth=2;
            elseif nargin==2
                linewidth=2;
                vectorsize=2;
            elseif nargin==3
                linewidth=2;
            end
            
            if isempty(obj.tr)
                obj=track_obj(obj);
                disp('Executed tracking first because tr was still empty')
            end
                        
            fname=[obj.foldername obj.filename '.tif'];
            I=imread(fname,frame);
                     
            figure
            hold on
            box on
            %show original image in grayscale
            colormap gray
            imagesc(I);
            sizeI=size(I);
            L=obj.NOP*vectorsize;
            %show vectorplot of droplet motion in red
            theta_t=obj.endvector(frame,2);
            x0=obj.tr_mass(frame,1);
            y0=obj.tr_mass(frame,2);
            x1=x0+L*cos(theta_t);
            x2=x0-L*cos(theta_t);
            y1=y0+L*sin(theta_t);
            y2=y0-L*sin(theta_t);
            x=linspace(x1,x2,10);
            y=y2+(y1-y2)/(x1-x2)*(x-x2);
            plot(x,y,'r-','LineWidth',linewidth)
            %Settings aspect ratios
            set(gca,'Ylim',[0 sizeI(1)],'Xlim', [0 sizeI(2)])
            %set(gca,'Ytick',[0:200:sizeI(2)],'Xtick',[0:300:sizeI(1)])
            pbaspect([sizeI(1)/sizeI(1) sizeI(1)/sizeI(2) 1]);
            %pbaspect([1 1 1])
            %axis([780 1280 300 800])
            set(gca,'YDir','normal');
            
        end
        function output=plot_DxDy_trend(obj,tmax)
            %plots Dx, Dy and Dbar (=0.5 Dtrans) as function of time
            %interval used for the measurement
           
            if isempty(obj.corrs)
                obj=calc_diff_tensor(obj);
            end
            if isempty(obj.MSD)
                obj=calc_MSD(obj);
            end
            if nargin<2 || tmax>obj.NOF
                tmax=obj.NOF;
            end
            
            dtmax=floor(tmax/10);
            Dxy=NaN(dtmax,3);
            %for each time interval
            for dt=1:dtmax
                %for each initial position
                Cxy=NaN(tmax-dt,2);
                for t=1:dt:tmax-dt
                    %the orientation angle in the intial positiston is
                    angle=-obj.endvector(t,2);
                    %gives a rotation matrix R
                    R=[cos(angle) -sin(angle) ; sin(angle) cos(angle)];
                    %get the trajectory for the object from time t onward
                    %and translate such that position at dt=0 is at origin
                    traj=obj.tr_mass(t+dt,1:2)-obj.tr_mass(t,1:2);
                    %rotate such that pos at time t is horizontal
                    traj_rot=R*traj';
                    %calculate the displacement in x and y
                    Cxy(t,1)=traj_rot(1,1)^2;
                    Cxy(t,2)=traj_rot(2,1)^2;
                end
                Dxy(dt,1)=dt;
                Dxy(dt,2)=(1/2)*nanmean(Cxy(:,1))/dt;
                Dxy(dt,3)=(1/2)*nanmean(Cxy(:,2))/dt;
            end
            
            temp1=store_cmass(obj);
            temp1=calc_MSD(temp1);
            MSD_t=temp1.MSD{1,2};
            Dt=[(1/4)*MSD_t(2:end,2)./MSD_t(2:end,1) MSD_t(2:end,1)];
            
            figure
            %hold on
            box on
            semilogx(Dt(:,2)/obj.framerate,Dt(:,1)*obj.scale^2*obj.framerate,'k-',Dxy(:,1)/obj.framerate,Dxy(:,3)*obj.scale^2*obj.framerate,'ro-',Dxy(:,1)/obj.framerate,Dxy(:,2)*obj.scale^2*obj.framerate,'bo-')
            %semilogx(Dxy(:,1)/obj.framerate,Dxy(:,3)*obj.scale^2*obj.framerate,'ro-')
            %semilogx(Dxy(:,1)/obj.framerate,Dxy(:,2)*obj.scale^2*obj.framerate,'bo-')
            %semilogx(Dxy(:,1)/obj.framerate,nanmean(Dxy(:,2:3),2)*obj.scale^2*obj.framerate,'go-')
            %xlabel('\Delta t (s)');
            %ylabel('D (\mum^2/s)');
            %hold off
            
%             figure
%             hold on
%             plot(Dt(:,2)/obj.framerate,Dt(:,1)*obj.scale^2*obj.framerate,'k-',Dxy(:,1)/obj.framerate,Dxy(:,3)*obj.scale^2*obj.framerate,'ro-',Dxy(:,1)/obj.framerate,Dxy(:,2)*obj.scale^2*obj.framerate,'bo-')
%             box on
%             hold off
            
            output=[Dxy Dt(:,1)];
        end
        function plot_DxDy_trend_angle(obj,anglerange,angleav)
            %WORKS ONLY FOR TRIMERS
            %plots Dx, Dy and Dbar (=0.5 Dtrans) as function of time
            %interval used for the measurement
            %
            %anglerange is width of accepted angles starting from average
            %angle that is accepted in radians. 2pi/3 means all angles.
            %
            %angleav is the average angle around which this is calculated
            %0 is stretched, 2 pi/3 is a closed triangle
            
           
            if isempty(obj.corrs)
                obj=calc_diff_tensor(obj);
            end
            if isempty(obj.MSD)
                obj=calc_MSD(obj);
            end
            if nargin<2
                anglerange=0.1;
            end
            if nargin<3
                angleav=0;
            end
            
            tmax=obj.NOF;
            dtmax=floor(tmax/10);
            Dxy=NaN(dtmax,3);
            %for each time interval
            for dt=1:dtmax
                %for each initial position
                Cxy=NaN(tmax-dt,2);
                for t=1:dt:tmax-dt
                    anglest=obj.angles{1,2};
                    angle_i=anglest(t,1);
                    %if the initial angle is in the right range
                    if abs(angle_i)>abs(angleav)-anglerange && abs(angle_i)<abs(angleav)+anglerange
                        %the orientation angle in the intial positiston is
                        angle=-obj.endvector(t,2);
                        %calculate the rotation matrix R
                        R=[cos(angle) -sin(angle) ; sin(angle) cos(angle)];
                        %get the trajectory for the object from time t onward
                        %and translate such that position at dt=0 is at origin
                        traj=obj.tr_mass(t+dt,1:2)-obj.tr_mass(t,1:2);
                        %rotate such that pos at time t is horizontal
                        traj_rot=R*traj';
                        %calculate the displacement in x and y
                        Cxy(t,1)=traj_rot(1,1)^2;
                        Cxy(t,2)=traj_rot(2,1)^2;
                    end
                end
                Dxy(dt,1)=dt;
                Dxy(dt,2)=(1/2)*nanmean(Cxy(:,1))/dt;
                Dxy(dt,3)=(1/2)*nanmean(Cxy(:,2))/dt;
            end
            
            temp1=store_cmass(obj);
            temp1=calc_MSD(temp1);
            MSD_t=temp1.MSD{1,2};
            Dt=[(1/4)*MSD_t(2:end,2)./MSD_t(2:end,1) MSD_t(2:end,1)];
            
            figure
            %hold on
            box on
            semilogx(Dt(:,2)/obj.framerate,Dt(:,1)*obj.scale^2*obj.framerate,'k-',Dxy(:,1)/obj.framerate,Dxy(:,3)*obj.scale^2*obj.framerate,'ro-',Dxy(:,1)/obj.framerate,Dxy(:,2)*obj.scale^2*obj.framerate,'bo-')
            %semilogx(Dxy(:,1)/obj.framerate,Dxy(:,3)*obj.scale^2*obj.framerate,'ro-')
            %semilogx(Dxy(:,1)/obj.framerate,Dxy(:,2)*obj.scale^2*obj.framerate,'bo-')
            %semilogx(Dxy(:,1)/obj.framerate,nanmean(Dxy(:,2:3),2)*obj.scale^2*obj.framerate,'go-')
            %xlabel('\Delta t (s)');
            %ylabel('D (\mum^2/s)');
            %hold off
            
%             figure
%             hold on
%             plot(Dt(:,2)/obj.framerate,Dt(:,1)*obj.scale^2*obj.framerate,'k-',Dxy(:,1)/obj.framerate,Dxy(:,3)*obj.scale^2*obj.framerate,'ro-',Dxy(:,1)/obj.framerate,Dxy(:,2)*obj.scale^2*obj.framerate,'bo-')
%             box on
%             hold off
            
            output=[Dxy Dt(:,1)];
        end
        function plot_scatter_correlations(obj,dt)
            %plots a scatter graph of the dx displacement vs dy
            %displacement for each timestep
            %Uses the reduced center of mass trajectory as input
            %dt is the timestep
            
            if nargin==1
                dt=1;
            end
            
            %check if required data are available
            if isempty(obj.tr_mass_red)
                obj=find_tr_mass_red(obj);
            end
            if isempty(obj.orients)
                obj=get_orients(obj);
            end
            if isempty(obj.angles)
                obj=find_angles(obj);
            end
            
            if obj.NOP>3
                disp('This code only works for a trimer right now')
            end
            
            %get displacements from center of mass trajectory
            xall=obj.tr_mass_red(:,1);
            x1=xall(1+dt:end);
            x2=xall(1:end-dt);
            xdisp=x1-x2;
            yall=obj.tr_mass_red(:,2);
            y1=yall(1+dt:end);
            y2=yall(1:end-dt);
            ydisp=y1-y2;
            
            %get angular displacement from orientation trajectory
            phi_all=obj.orients;
            phi1=obj.orients(1+dt:end,2);
            phi2=obj.orients(1:end-dt,2);
            phi_disp=phi1-phi2;
            
            %get internal angle displacement from angular trajectory
            alpha_all=obj.angles{2};
            alpha1=alpha_all(1+dt:end,1);
            alpha2=alpha_all(1:end-dt);
            alpha_disp=alpha1-alpha2;
            
            figure
            hold on
            box on
            plot(xdisp,ydisp,'k.');
            xlabel('\Deltax');
            ylabel('\Deltay');
            axis([-2.5 2.5 -2.5 2.5])
            pbaspect([1 1 1]);
            hold off
            
            figure
            hold on
            box on
            plot(xdisp,phi_disp,'k.');
            xlabel('\Deltax');
            ylabel('\Delta\theta');
            hold off
            
            figure
            hold on
            box on
            plot(ydisp,phi_disp,'k.');
            xlabel('\Deltay');
            ylabel('\Delta\theta');
            hold off
            
            figure
            hold on
            box on
            plot(xdisp,alpha_disp,'k.');
            xlabel('\Deltax');
            ylabel('\Delta\alpha');
            hold off
            
            figure
            hold on
            box on
            plot(ydisp,alpha_disp,'k.');
            xlabel('\Deltay');
            ylabel('\Delta\alpha');
            hold off
            
            figure
            hold on
            box on
            plot(phi_disp,alpha_disp,'k.');
            xlabel('\Delta\theta');
            ylabel('\Delta\alpha');
            hold off
            
            
        end
        function obj=crop_in_time(obj,tbegin,tend)
            %Returns the object with only data in between tbegin and tend
            %The only items it crops are obj.pos and obj.tr;
            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            if tbegin<1
                disp('Begin time cannot be smaller than 1.')
                tbegin=1;
            end
            
          tr_all=[];
          for p=1:obj.NOP
              tr_p=obj.tr(find(obj.tr(:,4)==p),1:4);
              nbegin=find(tr_p(:,3)==tbegin);
              nend=find(tr_p(:,3)==tend);
              tr_p=tr_p(nbegin:nend,:);
              tr_all=[tr_all; tr_p];
          end
          obj.tr=tr_all;
          obj.NOF=length(tr_p);
          obj=make_dist_matrix(obj);
          obj=find_neighbours(obj);
          obj=find_reduced_tracks(obj);
          obj=find_center_of_mass(obj);
          obj=find_angles(obj);
          
        end
        function output_temp=average_correlations(list_of_objects,imagesupress)
            %function takes in a list of objects of type single chain
            %organized as [obj_1 obj_2 obj_3... obj_n]
            %It outputs a new object of type single chain that contains the
            %average values of the individual objects.
            %all objects must have the same number of particles
            %all objects must have the same framerate and scalebar
            %imagesupress allows to supress the showing of images (0 is
            %show images, 1 is not show images). Default is 0 (show).
            
            if nargin==1
                imagesupress=0;
            end
            
            fname=inputname(1);
            
            %make sure the list is correctly formatted
            if length(list_of_objects)<2
                error('Not enough objects provided to average. Check input.');
            end
            
            %make sure the objects are of the right type
            object_1=list_of_objects(1);
            try NOP_1=object_1.NOP;
            catch
                error('Your list seems to not contain objects of type single chain. Please check input.');
            end
            framerate1=object_1.framerate;
            scale1=object_1.scale;
            
            %produce a new cell containing the correlation function of the
            %first particle
            all_corrs=object_1.corrs;
            all_diffs=object_1.diffs;
            all_angles=cell(NOP_1-2);
            for i=1:NOP_1+1
                for j=1:NOP_1+1
                    Cij=all_corrs{i,j}(:,2);
                    Dij=all_diffs{i,j};
                    if (i==1 || i==2) && (j==1 || j==2)
                        all_corrs{i,j}=Cij*scale1^2;
                        all_diffs{i,j}=Dij;
                    elseif i==1 || i==2 || j==1 || j==2
                        all_corrs{i,j}=Cij*scale1;
                        all_diffs{i,j}=Dij;
                    else
                        all_corrs{i,j}=Cij;
                        all_diffs{i,j}=Dij;
                    end
                end
            end
            for q=1:NOP_1-2
                all_angles{q,1}=q;
                all_angles{q,2}=object_1.angles{q,2};
            end
            %for all consecutive particles, append if framerate and
            %scalebar are okay
            for p=2:length(list_of_objects)
                %Test if correlations have been calculated for this object
                object_p=list_of_objects(p);
                NOP_p=object_p.NOP;
                %only proceed if object has same number of particles as
                %first object
                if NOP_p==NOP_1 && object_p.framerate==framerate1
                    
                    %fill the cell of correlations
                    corrs_p=object_p.corrs;
                    diffs_p=object_p.diffs;
                    for i=1:NOP_1+1
                        for j=1:NOP_1+1
                            minlength=min([length(all_corrs{i,j}(:,1)) length(corrs_p{i,j}(:,1))]);
                            Cij=corrs_p{i,j}(1:minlength,2);
                            Dij=diffs_p{i,j};
                            if (i==1 || i==2) && (j==1 || j==2)
                                all_corrs{i,j}=[all_corrs{i,j}(1:minlength,:) Cij*object_p.scale^2];    
                            elseif i==1 || i==2 || j==1 || j==2
                                all_corrs{i,j}=[all_corrs{i,j}(1:minlength,:) Cij*object_p.scale];
                            else
                                all_corrs{i,j}=[all_corrs{i,j}(1:minlength,:) Cij];
                            end
                            all_diffs{i,j}=[all_diffs{i,j} Dij];
                        end
                    end
                    for q=1:NOP_1-2
                        all_angles{q,2}=[all_angles{q,2} ; object_p.angles{q,2}];
                    end
                    
                else
                    disp(['Object ' num2str(p) ' in your list is ignored because it has ' num2str(object_p.NOP) ' particles and particle 1 has ' num2str(NOP_1) '.'])
                    disp('Alternatively the framerate might not be consistent over all your objects. Please check input.');
                end
            end
            
            %average correlations
            average_corrs=cell(NOP_1+1,NOP_1+1);
            std_corrs=cell(NOP_1+1,NOP_1+1);
            average_diffs=cell(NOP_1+1,NOP_1+1);
            average_diffs_2=cell(NOP_1+1,NOP_1+1);
            stds_diffs=cell(NOP_1+1,NOP_1+1);
            for i=1:NOP_1+1
                for j=1:NOP_1+1
                    Cij=[object_1.corrs{i,j}(1:length(all_corrs{i,j}(:,1)),1) mean(all_corrs{i,j},2)];
                    Cijstd=std(all_corrs{i,j},0,2);
                    average_corrs{i,j}=Cij;
                    std_corrs{i,j}=Cijstd;
                    Dij2=mean(all_diffs{i,j});
                    SDij=std(all_diffs{i,j});
                    Dij=(Cij(2,2)-Cij(1,2))/(Cij(2,1)-Cij(1,1));
                    average_diffs{i,j}=Dij*object_1.framerate;
                    average_diffs_2{i,j}=Dij2;
                    stds_diffs{i,j}=SDij;
                 end
            end
                        
            %Plot the average correlation functions
            f = figure(10);
            p = uipanel('Parent',f,'BorderType','none','BackgroundColor','white');
            p.Title = 'Diffusion tensor';
            p.TitlePosition = 'centertop';
            p.FontSize = 12;
            p.FontWeight = 'bold';
            for i=1:NOP_1+1
                for j=1:NOP_1+1
                    Cij=average_corrs{i,j};
                    Cijstd=std_corrs{i,j};
                    if imagesupress==0
                        hold on
                        box on
                        axi=subplot(NOP_1+1,NOP_1+1,(j-1)*(NOP_1+1)+i,'Parent',p);
                        errorbar(axi,Cij(:,1)/object_1.framerate,Cij(:,2),Cijstd(:,1),'k.');
                        title(['C' num2str(i) num2str(j)])
                        xlabel('dt (s)')
                        if (i==1 || i==2) && (j==1 || j==2)
                            ylabel('C(t) (\mum^2)');
                        elseif i==1 || i==2 || j==1 || j==2
                            ylabel('C(t) (\mum.rad)');
                        else
                            ylabel('C(t) (rad^2)');
                        end
                        hold off
                    end
                end
            end
            if imagesupress~=0
                close 10
            end
            
            %Create a new object with the averaged values
            parameters=get_parameters(object_1);

            output=single_chain(fname,parameters,'n');
            output.framerate=object_1.framerate;
            output.scale=object_1.scale;
            output.chainlength=object_1.chainlength;
            output.corrs=average_corrs;
            output.diffs=average_diffs;
            output.angles=all_angles;
            
            plot_angle_histogram(output);
            
            output_temp=cell(1,3);
            %single chain object with everything you might want
            output_temp{1}=output;
            %average diffusion coefficents calculated from individual diffs 
            output_temp{2}=average_diffs_2;
            %standard deviation in diffs
            output_temp{3}=stds_diffs;
            
            
            
                        
        end
        function plot_MSD_xy(obj)
            %plots the separate MSD functions in x direction, y direction
            %and the rotational average
            
            if isempty(obj.corrs)
                obj=calc_diff_tensor(obj);
            end
            
            Cxx=obj.corrs{1,1};
            Cyy=obj.corrs{2,2};
            Cav=[Cxx(:,1) 0.5*(Cxx(:,2)+Cyy(:,2))];
            
            figure
            hold on
            box on
            plot(Cxx(:,1)/obj.framerate,Cxx(:,2)*obj.scale^2,'b-');
            plot(Cyy(:,1)/obj.framerate,Cyy(:,2)*obj.scale^2,'r-');
            plot(Cav(:,1)/obj.framerate,Cav(:,2)*obj.scale^2,'k-');
            xlabel('\Delta t (s)');
            ylabel('MSD (\mum^2)');
            hold off
            
        end
    end
end