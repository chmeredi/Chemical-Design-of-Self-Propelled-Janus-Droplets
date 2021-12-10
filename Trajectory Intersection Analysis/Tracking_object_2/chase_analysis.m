classdef chase_analysis
    %Analysis code of movies of two types of droplets that chase each other
    %Designed to semi-automatically extract system parameters mobility,
    %activity and Peclet numbers 1 and 2.
    %
    %Contains a set of objects of type dropletpair<droplets<particles
    %The dropletpairs MUST consist of one encounter between two particles
    %only. This code bins and averages interactions between multiple encounters
    %to obtain statistically relevant numbers
    %
    %Functions
    %
    %Calculations
    %   - obj=chase_analysis(listofdropletpairs)         Initialization function
    %   - obj=add_dropletpair(obj,dropletpair)
    %   - obj=remove_dropletpair(obj,p)
    %   - output=find_minmaxdist(obj,n)
    %   - obj=average_u_rel(obj,binsize,n)
    %   - obj=average_u_sum(obj,binsize,n)
    %   - obj=average_ui(obj,binsize,n)
    %   - obj=extract_values(obj,lowlim,upperlim);
    %
    %Display functions
    %   - plot_u_rel_average(obj,plotlims)
    %   - plot_u_sum_average(obj,plotlims)
    %   - plot_ui_average(obj,plotlims)
    %   - plot_u_rel_rs(obj,plotlims)
    %   - plot_u_sum_rs(obj,plotlims)
    %   - plot_ui_rs(obj,plotlims)
    %
    %   - plot_u_rel_all(obj);
    %
    %
    %Properties
    %
    %   - scale
    %   - framerate
    %   - list_of_droplets
    %   - number_of_droplets
    %   - u_rel_average
    %   - u_sum_average
    %   - u1_average
    %   - u2_average
    %   - fitvals
    %
    %
    %Modification history
    %written on 19-07-2018 by Pepijn Moerman
    %
    
    properties
        scale
        framerate
        list_of_droplets
        number_of_droplets
        u_rel_average
        u_sum_average
        u1_average
        u2_average
        fitvals
    end
        methods
            function obj=chase_analysis(listofdropletpairs)
                %list of dropletpairs must be a 1xn list [obj1 obj2...]
                %of objects of type dropletpair
                %MAKE SURE YOU HAVE COMPLETED ANALYSIS ON INDIVIDUAL
                %ENCOUNTERS BEFORE ADDING TO THIS CODE
                
                obj.list_of_droplets={};
                
                for p=1:length(listofdropletpairs)
                    droplet_p=listofdropletpairs(p);
                    %droplet_p=analyze(droplet_p);
                    obj.list_of_droplets{p}=droplet_p;
                end
                obj.number_of_droplets=length(listofdropletpairs);
                obj.scale=droplet_p.scale;
                obj.framerate=droplet_p.framerate;
                
            end
            function obj=change_scale(obj,scalenew)
                %change the scale to new value
                obj.scale=scalenew;
            end
            function obj=change_framerate(obj,frameratenew)
                %change the scale to new value
                obj.framerate=frameratenew;
            end
            function obj=add_dropletpair(obj,dropletpair)
                %adds an instance of dropletpair to the object
                
                listtemp=obj.list_of_droplets;
                listtemp{end+1}=dropletpair;
                obj.list_of_droplets=listtemp;
                obj.number_of_droplets=obj.number_of_droplets+1;
                
            end
            function obj=shift_u_sum(obj,amount)
                %If u_rel does not decay to zero at large separation, this
                %function corrects by shifting the whole function a chosen
                %amount
                %amount is added to u_rel
                %give amount in micrometer/second
                
                obj.u_sum_average(:,2)=obj.u_sum_average(:,2)+amount;
                     
            end
            function obj=shift_u_rel(obj,amount)
                %If u_rel does not decay to zero at large separation, this
                %function corrects by shifting the whole function a chosen
                %amount
                %amount is added to u_rel
                %give amount in micrometer/second
                
                obj.u_rel_average(:,2)=obj.u_rel_average(:,2)+amount;
                     
            end
            function obj=shift_u1(obj,amount)
                %If u_rel does not decay to zero at large separation, this
                %function corrects by shifting the whole function a chosen
                %amount
                %amount is added to u_rel
                %give amount in micrometer/second
                
                obj.u1_average(:,2)=obj.u1_average(:,2)+amount;
                     
            end
            function obj=shift_u2(obj,amount)
                %If u_rel does not decay to zero at large separation, this
                %function corrects by shifting the whole function a chosen
                %amount
                %amount is added to u_rel
                %give amount in micrometer/second
                
                obj.u2_average(:,2)=obj.u2_average(:,2)+amount;
                     
            end
            function obj=remove_dropletpair(obj,p)
                %removes particle with index pfrom the list of dropletpairs
                                
                if p<=obj.number_of_droplets
                    listnew={};
                    count=0;
                    listtemp=obj.list_of_droplets;
                    for q=1:obj.number_of_droplets
                        if q~=p
                            count=count+1;
                            listnew{count}=listtemp{q};
                        end
                    end
                else
                    disp(['There is no droplet with index ' num2str(p) '.'])
                end
                obj.list_of_droplets=listnew;
                obj.number_of_droplets=obj.number_of_droplets-1;
                
            end
            function output=find_minmaxdist(obj,n)
                %Finds the minimum and maximum distance between particles with at least n droplets
                %default n=2; n can also not be smaller than 2. Then
                %averages and standard deviations are a problem.
                
                %check that there are at least 2 droplets present
                if obj.number_of_droplets<2
                    error('This function only works if there are at least 2 instances of dropletpair2 saved in the object.');
                elseif obj.number_of_droplets==2
                    n=1;
                end
                
                %give default input
                if nargin<2
                    n=2;
                end
                if n<2
                    n=1;
                end
                
                %give initial values min and max dist
                mindist=NaN(obj.number_of_droplets,1);
                maxdist=NaN(obj.number_of_droplets,1);
                
                %loop through particle list and find smallest and largest r
                for p=1:obj.number_of_droplets
                    %unwrap instance of dropletpair2
                    droplet_p=obj.list_of_droplets{p};
                    distance_p=droplet_p.distance(:,2)*obj.scale;
                    maxdist(p,1)=max(distance_p(:,1));
                    mindist(p,1)=min(distance_p(:,1));
                end
                
                %sort the list of minimal distances
                mindist(:,1)=sort(mindist(:,1));
                maxdist(:,1)=sort(maxdist(:,1));
                
                %take the nth element of the list. That will be a distance
                %that at least n droplets have data for
                mindistall=mindist(n,1);
                maxdistall=maxdist(end-n,1);
                
                %package output
                output=[mindistall maxdistall];
                
            end
            function obj=average_u_rel(obj,binsize,n)
                % Bins the v_rel (interactionstrength) for all instances of
                % dropletpair in the object using the same binsize
                % "binssize"
                % binsize is given in micrometers. Data are saved in
                % micrometers/second
                % n is the number of datasets you want to have left at
                % least to average over
                
                %make sure enought arguments are provided
                if nargin<3
                    n=2;
                end
                if nargin<2
                    binsize=1;
                end
                                
                %get minimum and maximum distance
                distlims=find_minmaxdist(obj,n);
                
                %makes sure binsize is a integer
                binsize=round(binsize);
                if binsize<1
                    binsize=1;
                elseif binsize>(distlims(2)-distlims(1))*obj.scale
                    error('Binsize too large. Please pick a smaller value.');
                end
                
                %add all data on u_rel to big list
                u_rel_all=[];
                for p=1:obj.number_of_droplets
                    droplet_p=obj.list_of_droplets{p};
                    u_rel_p=droplet_p.u_rel(:,1:2);
                    u_rel_p(:,1)=u_rel_p(:,1)*obj.scale;
                    u_rel_p(:,2)=u_rel_p(:,2)*obj.scale*droplet_p.framerate;
                    u_rel_all=[u_rel_all; u_rel_p];
                end
                
                %get u_rel_all real units (um and um/s);
                u_rel_all(:,1)=u_rel_all(:,1);%*obj.scale;
                u_rel_all(:,2)=u_rel_all(:,2);%*obj.scale*obj.framerate;
                
                %get the minimum and maximum distance in whole values in
                %micrometers
                dmin=ceil(distlims(1));%;*obj.scale);
                dmax=floor(distlims(2));%;*obj.scale);
                
                %create an empty matrix for average speed of the right size
                nbins=round((dmax-dmin)/binsize+1);
                u_rel_av=NaN(nbins,3);
                
                %loop through all u_rel data and average per bin
                for n=1:nbins
                    
                    %get the bin limits for bin n
                    binlimlow=dmin+(n-1)*binsize;
                    binlimhigh=dmin+(n)*binsize;
                    
                    %save the distance in the first column
                    u_rel_av(n,1)=binlimlow;
                    
                    %get the average and stdev in speed and save in second and third column;
                    u_rel_n=u_rel_all(find(u_rel_all(:,1)>=binlimlow & u_rel_all(:,1)<binlimhigh),2);
                    u_rel_n_av=nanmean(u_rel_n(:,1));
                    u_rel_n_std=nanstd(u_rel_n(:,1));
                    u_rel_av(n,2)=u_rel_n_av;
                    u_rel_av(n,3)=u_rel_n_std;
                    
                end
                
                obj.u_rel_average=u_rel_av;
                
            end
            function obj=average_u_sum(obj,binsize,n)
                % Bins the v_sum chase speed for all instances of
                % dropletpair in the object using the same binsize
                % "binssize"
                % binsize is given in micrometers. Data are saved in
                % micrometers/second
                
                %make sure enought arguments are provided
                if nargin<3
                    n=2;
                end
                if nargin<2
                    binsize=1;
                end
                                                
                %get minimum and maximum distance
                distlims=find_minmaxdist(obj,n);
                
                %makes sure binsize is a integer
                binsize=round(binsize);
                if binsize<1
                    binsize=1;
                elseif binsize>(distlims(2)-distlims(1))*obj.scale
                    error('Binsize too large. Please pick a smaller value.');
                end
                
                %add all data on u_rel to big list
                u_sum_all=[];
                for p=1:obj.number_of_droplets
                    droplet_p=obj.list_of_droplets{p};
                    u_sum_p=droplet_p.u_sum(:,1:2);
                    u_sum_p(:,1)=u_sum_p(:,1)*obj.scale;
                    u_sum_p(:,2)=u_sum_p(:,2)*obj.scale*droplet_p.framerate;
                    u_sum_all=[u_sum_all; u_sum_p];
                end
                
                %get u_rel_all real units (um and um/s);
                u_sum_all(:,1)=u_sum_all(:,1);%*obj.scale;
                u_sum_all(:,2)=u_sum_all(:,2);%*obj.scale*obj.framerate;
                
                %get the minimum and maximum distance in whole values in
                %micrometers
                dmin=ceil(distlims(1));%;*obj.scale);
                dmax=floor(distlims(2));%;*obj.scale);
                
                %create an empty matrix for average speed of the right size
                nbins=round((dmax-dmin)/binsize+1);
                u_sum_av=NaN(nbins,3);
                
                nbins
                %loop through all u_rel data and average per bin
                for n=1:nbins
                    
                    %get the bin limits for bin n
                    binlimlow=dmin+(n-1)*binsize;
                    binlimhigh=dmin+(n)*binsize;
                    
                    %save the distance in the first column
                    u_sum_av(n,1)=binlimlow;
                    
                    %get the average and stdev in speed and save in second and third column;
                    u_sum_n=u_sum_all(find(u_sum_all(:,1)>=binlimlow & u_sum_all(:,1)<binlimhigh),2);
                    u_sum_n_av=nanmean(u_sum_n(:,1));
                    u_sum_n_std=nanstd(u_sum_n(:,1));
                    u_sum_av(n,2)=u_sum_n_av;
                    u_sum_av(n,3)=u_sum_n_std;
                    
                end
                
                obj.u_sum_average=u_sum_av;
                
            end
            function obj=average_ui(obj,binsize,n)
                % Bins the v_sum chase speed for all instances of
                % dropletpair in the object using the same binsize
                % "binssize"
                % binsize is given in micrometers. Data are saved in
                % micrometers/second
                
                %make sure enought arguments are provided
                if nargin<3
                    n=2;
                end
                if nargin<2
                    binsize=1;
                end
                                
                %get minimum and maximum distance
                distlims=find_minmaxdist(obj,n);
                
                %makes sure binsize is a integer
                binsize=round(binsize);
                if binsize<1
                    binsize=1;
                elseif binsize>(distlims(2)-distlims(1))*obj.scale
                    error('Binsize too large. Please pick a smaller value.');
                end
                
                %add all data on u1 and u2 to big list
                u1_all=[];
                u2_all=[];
                for p=1:obj.number_of_droplets
                    droplet_p=obj.list_of_droplets{p};
                    u1_p=droplet_p.u1(:,1:2);
                    u2_p=droplet_p.u2(:,1:2);
                    u1_p(:,1)=u1_p(:,1)*obj.scale;
                    u1_p(:,2)=u1_p(:,2)*obj.scale*droplet_p.framerate;
                    u2_p(:,1)=u2_p(:,1)*obj.scale;
                    u2_p(:,2)=u2_p(:,2)*obj.scale*droplet_p.framerate;
                    u1_all=[u1_all; u1_p];
                    u2_all=[u2_all; u2_p];
                end
                
                %get u1_all and u2_all real units (um and um/s);
                u1_all(:,1)=u1_all(:,1);%*obj.scale;
                u1_all(:,2)=u1_all(:,2);%*obj.scale*obj.framerate;
                u2_all(:,1)=u2_all(:,1);%*obj.scale;
                u2_all(:,2)=u2_all(:,2);%*obj.scale*obj.framerate;
                
                %get the minimum and maximum distance in whole values in
                %micrometers
                dmin=ceil(distlims(1));%*obj.scale);
                dmax=floor(distlims(2));%*obj.scale);
                
                %create empty matrices for average speed of the right size
                nbins=round((dmax-dmin)/binsize+1);
                u1_av=NaN(nbins,3);
                u2_av=NaN(nbins,3);
                
                %loop through all u_rel data and average per bin
                for n=1:nbins
                    
                    %get the bin limits for bin n
                    binlimlow=dmin+(n-1)*binsize;
                    binlimhigh=dmin+(n)*binsize;
                    
                    %save the distance in the first column
                    u1_av(n,1)=binlimlow;
                    u2_av(n,1)=binlimlow;
                    
                    %get the average and stdev in speed u1 and save in second and third column;
                    u1_n=u1_all(find(u1_all(:,1)>=binlimlow & u1_all(:,1)<binlimhigh),2);
                    u1_n_av=nanmean(u1_n(:,1));
                    u1_n_std=nanstd(u1_n(:,1));
                    u1_av(n,2)=u1_n_av;
                    u1_av(n,3)=u1_n_std;
                    
                    %get the average and stdev in speed u2 and save in second and third column;
                    u2_n=u2_all(find(u2_all(:,1)>=binlimlow & u2_all(:,1)<binlimhigh),2);
                    u2_n_av=nanmean(u2_n(:,1));
                    u2_n_std=nanstd(u2_n(:,1));
                    u2_av(n,2)=u2_n_av;
                    u2_av(n,3)=u2_n_std;
                    
                end
                
                %save average to object
                obj.u1_average=u1_av;
                obj.u2_average=u2_av;
                
            end
            function obj=extract_values(obj,lowerlim,upperlim)
                %Extracts the interaction strength of two droplets by
                %fitting the U1 vs 1/r2 and U2 vs 1/r2
                %the prefactor for U1 is called C12 and U2 is called C21.
                
                %check if input exists
                if isempty(obj.u1_average)
                    obj=average_ui(obj);
                end
                
                %obtain input matrices
                u1_av=obj.u1_average;
                u2_av=obj.u2_average;
                
                %make sure enough arguments are provided
                if nargin<3
                    upperlim=max(u1_av(:,1));
                end
                if nargin<2
                    lowerlim=min(u1_av(:,1));
                end
                
                %find right range so that the range falls within both the
                %user provided limits and the measured distance range
                rmin=max([lowerlim min(u1_av(:,1))]);
                rmax=min([upperlim max(u1_av(:,1))]);
                
                %clip data to right range
                u1_clip=u1_av(u1_av(:,1)>=rmin & u1_av(:,1)<=rmax,:);
                u2_clip=u2_av(u1_av(:,1)>=rmin & u2_av(:,1)<=rmax,:);
                u1_clip(isnan(u1_clip(:,2)),:)=[];
                u2_clip(isnan(u2_clip(:,2)),:)=[];
                %u1_clip=u1_av(find(u1_av(:,1)==rmin):find(u1_av(:,1)==rmax),:);
                %u2_clip=u2_av(find(u2_av(:,1)==rmin):find(u2_av(:,1)==rmax),:);
                                
                %fit data
                [C12,S12]=polyfit(u1_clip(:,1).^(-2),u1_clip(:,2),1);
                [C21,S21]=polyfit(u2_clip(:,1).^(-2),u2_clip(:,2),1);
               
                %plot data u1
                x=u1_clip(:,1).^(-2);
                y=u1_clip(:,2);
                error=u1_clip(:,3);
                x1=linspace(rmin-3,rmax+3,100).^(-2);
                y1=polyval(C12,x1);
                figure
                hold on
                box on
                title('u1')
                errorbar(x,y,error,'ksq')
                plot(x1,y1,'r-');                
                xlabel('r^{-2} (\mum^{-2})');
                ylabel('v (\mum/s)');
                axis([0.1*10^-4 20.0*10^-4 -20 20])
                yticks([-20 -10 0 10 20])
                hold off
                
                %plot data u2
                x=u2_clip(:,1).^(-2);
                y=u2_clip(:,2);
                error=u2_clip(:,3);
                x1=linspace(rmin-3,rmax+3,100).^(-2);
                y1=polyval(C21,x1);
                figure
                hold on
                box on
                title('u2')
                errorbar(x,y,error,'ksq')
                plot(x1,y1,'r-');                
                xlabel('r^{-2} (\mum^{-2})');
                ylabel('v (\mum/s)');
                axis([0.1*10^-4 20.0*10^-4 -20 20])
                yticks([-20 -10 0 10 20])
                hold off
                
                %Calculating the covariance matrix from fit as Rinv*Rinv')*normr^2/df;
                try p12=inv(S12.R)*inv(S12.R)'*S12.normr^2/S12.df;
                catch
                    p12=[NaN,NaN];
                    disp('Could not determine fit accuracy for this limit')
                end
                try p21=inv(S21.R)*inv(S21.R)'*S21.normr^2/S21.df;
                catch
                    p21=[NaN,NaN];
                    disp('Could not determine fit accuracy for this limit')
                end
                
                %std fit val C12 = p12(1,2);
                %save the values as C12 and there stdevs and C21 (where C12 is slope of speed of U1)
                obj.fitvals=[C12(1) p12(1,2); C21(1) p21(1,2)];
                 
            end
            
            %-------------------- DISPLAY FUNCTIONS ----------------------%
            function plot_u_relsum_bb(obj,plotlims)
                %plots the average u_rel vs r
                %plotlims are the plot boundaries [xmin xmax ymin ymax];
                
                %check if input exists
                if isempty(obj.u_rel_average)
                    obj=average_u_rel(obj);
                end
                
                %get x y and errorbar values
                x=obj.u_rel_average(:,1)-obj.u_rel_average(1,1);
                y=obj.u_rel_average(:,2);
                error=obj.u_rel_average(:,3);
                
                %plot the figure
                figure
                hold on
                box on
                errorbar(x,y,error,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r (\mum)');
                ylabel('v (\mum/s)');
                hold off
                
                %plots the average u_sum vs r
                %plotlims are the plot boundaries [xmin xmax ymin ymax];
                
                %check if input exists
                if isempty(obj.u_sum_average)
                    obj=average_u_sum(obj);
                end
                
                %get x y and errorbar values
                x=obj.u_sum_average(:,1)-obj.u_sum_average(1,1);
                y=obj.u_sum_average(:,2);
                error=obj.u_sum_average(:,3);
                
                %plot the figure
                figure
                hold on
                box on
                errorbar(x,y,error,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r (\mum)');
                ylabel('v (\mum/s)');
                hold off
            end
            function plot_u_both_bb(obj)
                %plots the average u_rel & u sum vs r
                %plotlims are the plot boundaries [xmin xmax ymin ymax];
                
                %check if input exists
                if isempty(obj.u_rel_average)
                    obj=average_u_rel(obj);
                end
                
                %get x y and errorbar values for urel
                x1=obj.u_rel_average(:,1)-obj.u_rel_average(1,1);
                y1=obj.u_rel_average(:,2);
                error1=obj.u_rel_average(:,3);
                
                %get x and y and errorbar values for usum
                x2=obj.u_sum_average(:,1)-obj.u_sum_average(1,1);
                y2=obj.u_sum_average(:,2);
                error2=obj.u_sum_average(:,3);
                
                %plot the figure
                figure
                hold on
                box on
                errorbar(x1,y1,error1,'kd')
                errorbar(x2,y2,error2,'rd')
                axis([-5 80 -10 10])
                yticks([-10 -5 0 5 10]);
                xticks([0 20 40 60 80]);
                xlabel('r (\mum)');
                ylabel('v (\mum/s)');
                hold off
                
                
            end
            function plot_u_rel_average(obj,plotlims)
                %plots the average u_rel vs r
                %plotlims are the plot boundaries [xmin xmax ymin ymax];
                
                %check if input exists
                if isempty(obj.u_rel_average)
                    obj=average_u_rel(obj);
                end
                
                %get x y and errorbar values
                x=obj.u_rel_average(:,1);
                y=obj.u_rel_average(:,2);
                error=obj.u_rel_average(:,3);
                
                %plot the figure
                figure
                hold on
                box on
                errorbar(x,y,error,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r (\mum)');
                ylabel('v (\mum/s)');
                hold off
            end
            function plot_u_sum_average(obj,plotlims)
                %plots the average u_sum vs r
                %plotlims are the plot boundaries [xmin xmax ymin ymax];
                
                %check if input exists
                if isempty(obj.u_sum_average)
                    obj=average_u_sum(obj);
                end
                
                %get x y and errorbar values
                x=obj.u_sum_average(:,1);
                y=obj.u_sum_average(:,2);
                error=obj.u_sum_average(:,3);
                
                %plot the figure
                figure
                hold on
                box on
                errorbar(x,y,error,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r (\mum)');
                ylabel('v (\mum/s)');
                hold off
            end
            function plot_ui_average(obj,plotlims)
                %plots the average u1 and u2 vs r
                %plotlims are the plot boundaries [xmin xmax ymin ymax];
                
                %check if input exists
                if isempty(obj.u1_average)
                    obj=average_ui(obj);
                end
                
                %get x y and errorbar values
                x=obj.u1_average(:,1);
                y=obj.u1_average(:,2);
                error=obj.u1_average(:,3);
                
                %plot the figure
                figure
                hold on
                box on
                title('U1')
                errorbar(x,y,error,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r (\mum)');
                ylabel('v (\mum/s)');
                hold off
                
                %get x y and errorbar values
                x=obj.u2_average(:,1);
                y=obj.u2_average(:,2);
                error=obj.u2_average(:,3);
                
                %plot the figure
                figure
                hold on
                box on
                title('U2')
                errorbar(x,y,error,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r (\mum)');
                ylabel('v (\mum/s)');
                hold off
            end
            function plot_u_sum_rs(obj,plotlims)
                %plots the average u_sum vs 1/r2
                %plotlims are the plot boundaries [xmin xmax ymin ymax];
                
                %check if input exists
                if isempty(obj.u_sum_average)
                    obj=average_u_sum(obj);
                end
                
                %get x y and errorbar values
                x=(obj.u_sum_average(:,1)).^(-2);
                y=obj.u_sum_average(:,2);
                error=obj.u_sum_average(:,3);
                
                %plot the figure
                figure
                hold on
                box on
                errorbar(x,y,error,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r^{-2} (\mum^{-2})');
                ylabel('v (\mum/s)');
                hold off
            end
            function plot_u_rel_rs(obj,plotlims)
                %plots the average u_rel vs 1/r2
                %plotlims are the plot boundaries [xmin xmax ymin ymax];
                
                %check if input exists
                if isempty(obj.u_rel_average)
                    obj=average_u_rel(obj);
                end
                
                %get x y and errorbar values
                x=(obj.u_rel_average(:,1)).^(-2);
                y=obj.u_rel_average(:,2);
                error=obj.u_rel_average(:,3);
                
                %plot the figure
                figure
                hold on
                box on
                errorbar(x,y,error,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r^{-2} (\mum^{-2})');
                ylabel('v (\mum/s)');
                hold off
            end
            function plot_ui_rs(obj,plotlims)
                %plots the average u_rel vs 1/r2
                %plotlims are the plot boundaries [xmin xmax ymin ymax];
                
                %check if input exists
                if isempty(obj.u1_average)
                    obj=average_ui(obj);
                end
                
                %get x y and errorbar values for u1
                x=(obj.u1_average(:,1)).^(-2);
                y=obj.u1_average(:,2);
                error=obj.u1_average(:,3);
                
                %plot the figure
                figure
                hold on
                title('u1')
                box on
                errorbar(x,y,error,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r^{-2} (\mum^{-2})');
                ylabel('v (\mum/s)');
                hold off
                
                %get x y and errorbar values for u2
                x=(obj.u2_average(:,1)).^(-2);
                y=obj.u2_average(:,2);
                error=obj.u2_average(:,3);
                
                %plot the figure
                figure
                hold on
                title('u2')
                box on
                errorbar(x,y,error,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r^{-2} (\mum^{-2})');
                ylabel('v (\mum/s)');
                hold off
            end
            function plot_u_rel_all(obj,plotlims)
                %plots the average u_sum vs 1/r2
                %plotlims are the plot boundaries [xmin xmax ymin ymax];
                
                u_rel_all=[];
                for p=1:obj.number_of_droplets
                    droplet_p=obj.list_of_droplets{p};
                    u_rel_all=[u_rel_all ; droplet_p.u_rel];
                end
                
                %get x y and errorbar values
                x=u_rel_all(:,1);
                y=u_rel_all(:,2);
                
                %plot the figure
                figure
                hold on
                box on
                plot(x,y,'ksq')
                if nargin>1
                    axis(plotlims)
                end
                xlabel('r (\mum)');
                ylabel('v (\mum/s)');
                hold off
            end
       end%end methods
end%end class