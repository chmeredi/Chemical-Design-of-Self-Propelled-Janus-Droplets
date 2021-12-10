classdef chase_conditions
    %Class designed to plot average interaction strength of chasing
    %droplets for different conditions in one plot
    %
    %Functions
    %
    %Properties
    %
    %
    %
    %
    %
    %
    properties
        listofconditions
        numberofconditions
        integralvalues
    end
    methods
        function obj=chase_conditions(listofchase_analysis)
            %initialization function
            %MAKE SURE YOU HAVE FINISHED ANALYSIS OF INDIVIDUAL CONDITIONS
            %BEFORE CREATING THIS OBJECT.
            %this object is for plotting purposes ONLY.
            
            obj.listofconditions={};
            
            for p=1:length(listofchase_analysis)
                droplet_p=listofchase_analysis(p);
                obj.listofconditions{p}=droplet_p;
            end
            
            obj.numberofconditions=p;
            
        end
        function plot_ui_loglog(obj)
        end
        function plot_ui_all(obj)
            %plots the speed of the chaser and chasee as function of
            %interparticle distance for all conditions
            
            %get colormap
            cmap=winter(obj.numberofconditions);
            
            %make figure u1
            figure
            hold on
            box on
            title('U1')
            for p=1:obj.numberofconditions
                %Get values to plot
                condip=obj.listofconditions{p};
                u1_p=condip.u1_average;
                x=u1_p(:,1);
                y=u1_p(:,2);
                error=u1_p(:,3);
                %plot values
                errorbar(x,y,error,'sq','Color',cmap(p,:))
            end
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            hold off
            
            %make figure u2
            figure
            hold on
            box on
            title('U2')
            for p=1:obj.numberofconditions
                %Get values to plot
                condip=obj.listofconditions{p};
                u2_p=condip.u2_average;
                x=u2_p(:,1);
                y=u2_p(:,2);
                error=u2_p(:,3);
                %plot values
                errorbar(x,y,error,'sq','Color',cmap(p,:))
            end
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            hold off
        end
        function plot_fit_lin_bb(obj)
            %plot the approach and escape speed of droplets on a lin lin
            %scale U vs r with fits of shape U(r) = A r^-2 + B
            %plots r as boundary-boundary distance so r=0 is contact
            
            cmap=winter(obj.numberofconditions);
            
            %plot figure 1
            figure
            hold on
            box on
            title('u1')
            for p=1:obj.numberofconditions
                
                condip=obj.listofconditions{p};
                
                %obtain input matrices
                u1_av=condip.u1_average;
                                               
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
                u1_clip=u1_av(find(u1_av(:,1)==rmin):find(u1_av(:,1)==rmax),:);
                u1_clip(isnan(u1_clip(:,2)),:)=[];
                
                %fit data
                [C12,S12]=polyfit(u1_clip(:,1).^(-2),u1_clip(:,2),1);
                                
                %plot data u1
                x=u1_clip(:,1)-u1_clip(1,1);
                y=u1_clip(:,2);
                error=u1_clip(:,3);
                x1=linspace(rmin-3,rmax+3,100);
                y1=C12(1)./x1.^2;
                y1=y1(:)-y1(100);
                errorbar(x,y,error,'sq','Color',cmap(p,:))
                plot(x1-rmin,y1,'-','Color',cmap(p,:));
                
                
            end
            xlabel('r (um)');
            ylabel('v (um/s)');
            axis([0 80 -20 20])
            yticks([-20 -10 0 10 20])
            xticks([0 20 40 60 80])
            hold off
            
            %plot figure 2
            figure
            hold on
            box on
            title('u2')
            for p=1:obj.numberofconditions
                
                condip=obj.listofconditions{p};
                
                %obtain input matrices
                u2_av=condip.u2_average;
                u1_av=condip.u1_average;
                
                %make sure enough arguments are provided
                if nargin<3
                    upperlim=max(u2_av(:,1));
                end
                if nargin<2
                    lowerlim=min(u2_av(:,1));
                end
                
                %find right range so that the range falls within both the
                %user provided limits and the measured distance range
                rmin=max([lowerlim min(u2_av(:,1))]);
                rmax=min([upperlim max(u2_av(:,1))]);
                
                %clip data to right range
                u2_clip=u2_av(find(u2_av(:,1)==rmin):find(u2_av(:,1)==rmax),:);
                u2_clip(isnan(u2_clip(:,2)),:)=[];
                
                %fit data
                [C21,S21]=polyfit(u2_clip(:,1).^(-2),u2_clip(:,2),1);
                             
                %plot data u2
                x=u2_clip(:,1)-u2_clip(1,1);
                y=u2_clip(:,2);
                error=u2_clip(:,3);
                x1=linspace(rmin-3,rmax+3,100);
                y1=C21(1)./x1.^2;
                y1=y1(:)-y1(100);
                errorbar(x,y,error,'sq','Color',cmap(p,:))
                plot(x1-rmin,y1,'-','Color',cmap(p,:));
                
                
            end
            xlabel('r^ (um)');
            ylabel('v (um/s)');
            axis([0 80 -20 20])
            yticks([-20 -10 0 10 20])
            xticks([0 20 40 60 80])
            hold off
        end
        function plot_fit_lin_both_bb(obj,lowerlim,upperlim)
            %plot the approach and escape speed of droplets on a lin lin
            %scale U vs r with fits of shape U(r) = A r^-2 + B
            %plots r as boundary-boundary distance so r=0 is contact
            
            cmap=winter(obj.numberofconditions);
            
            %plot figure 1
            figure
            hold on
            box on
            for p=1:obj.numberofconditions
               
                condip=obj.listofconditions{p};
                
                %obtain input matrices
                u1_av=condip.u1_average;
                u2_av=condip.u2_average;
                                               
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
                u1_clip=u1_av(find(u1_av(:,1)==rmin):find(u1_av(:,1)==rmax),:);
                u1_clip(isnan(u1_clip(:,2)),:)=[];
                
                %clip data to right range
                u2_clip=u2_av(find(u2_av(:,1)==rmin):find(u2_av(:,1)==rmax),:);
                u2_clip(isnan(u2_clip(:,2)),:)=[];
                
                %fit data
                [C12,S12]=polyfit(u1_clip(:,1).^(-2),u1_clip(:,2),1);
                
                %fit data
                [C21,S21]=polyfit(u2_clip(:,1).^(-2),u2_clip(:,2),1);
                                
                %plot data u1
                x=u1_av(:,1)-u1_av(1,1);
                y=u1_av(:,2);
                error=u1_av(:,3);
                x1=linspace(rmin,rmax,100);
                y1=C12(1)./x1.^2;
                y1=y1(:)-y1(end);
                errorbar(x,y,error,'sq','Color',cmap(p,:))
                plot(x1-rmin,y1,'-','Color',cmap(p,:));
                
                                            
                %plot data u2
                x2=u2_av(:,1)-u2_av(1,1);
                y2=u2_av(:,2);
                error=u2_av(:,3);
                x3=linspace(rmin,rmax,100);
                y3=C21(1)./x3.^2;
                y3=y3(:)-y3(end);
                errorbar(x2,y2,error,'o','Color',cmap(p,:))
                plot(x3-rmin,y3,'-','Color',cmap(p,:));
                
            end
            xlabel('r (um)');
            ylabel('v (um/s)');
            axis([-5 60 -18 18])
            yticks([-15 -10 -5 0 5 10 15])
            xticks([0 20 40 60])
            hold off
            
           
        end
        function plot_maxval(obj,listofxvals)
            %plots the maximum valuues of the ui vs r curves
            %as a function of the control parameter
            
            u1_maxvals=NaN(length(listofxvals),3);
            u2_maxvals=NaN(length(listofxvals),3);
            u1_maxvals(:,1)=listofxvals';
            u2_maxvals(:,1)=listofxvals';
            
            figure
            hold on
            box on
            for p=1:obj.numberofconditions
                condip=obj.listofconditions{p};
                u1_p=condip.u1_average;
                u2_p=condip.u2_average;
                
                u1_high=max(u1_p(:,2));
                u1_low=min(u1_p(:,2));
                if abs(u1_low)>abs(u1_high)
                    u1_max=u1_low;
                else
                    u1_max=u1_high;
                end
                
                u2_high=max(u2_p(:,2));
                u2_low=min(u2_p(:,2));
                if abs(u2_low)>abs(u2_high)
                    u2_max=u2_low;
                else
                    u2_max=u2_high;
                end
                
                %u1_max=max(u1_p(:,2));
                %u2_max=max(u2_p(:,2));
                u1_error=u1_p(find(u1_p(:,2)==u1_max),3);
                u2_error=u2_p(find(u2_p(:,2)==u2_max),3);
                u1_maxvals(p,2:3)=[u1_max u1_error];
                u2_maxvals(p,2:3)=[u2_max u2_error];
            end
            errorbar(u1_maxvals(:,1),u1_maxvals(:,2),u1_maxvals(:,3),'rd','MarkerSize',12);
            errorbar(u2_maxvals(:,1),u2_maxvals(:,2),u2_maxvals(:,3),'bd','MarkerSize',12);
            xlabel('xvalues (scale)');
            ylabel('v_{max} (um/s)');
            axis([0 6 -20 20])
            xticks=([6 8 10 12 14 16]);
            yticks([-20 -10 0 10 20]);
            hold off
        end
        function obj=calc_integral_values(obj,listofxvals)
            %calculates and plots the integral under the u1 and u2 curves
            %of the ui vs r curves as a function of the control parameter
            %The interpretation of the values is the dissipated energy
            
            u1_intvals=NaN(length(listofxvals),3);
            u2_intvals=NaN(length(listofxvals),3);
            u1_intvals(:,1)=listofxvals';
            u2_intvals(:,1)=listofxvals';
            
            figure
            hold on
            box on
            for p=1:obj.numberofconditions
                condip=obj.listofconditions{p};
                u1_p=condip.u1_average;
                u2_p=condip.u2_average;
                
                
                u1_int=nansum(u1_p(:,2));
                u2_int=nansum(u2_p(:,2));
                
                u1_error=sqrt(nansum(u1_p(:,3).^2));
                u2_error=sqrt(nansum(u2_p(:,3).^2));
                %u1_error=nansum(u1_p(:,3));
                %u2_error=nansum(u2_p(:,3));
                u1_intvals(p,2:3)=[u1_int u1_error];
                u2_intvals(p,2:3)=[u2_int u2_error];
            end
            errorbar(u1_intvals(:,1),u1_intvals(:,2),u1_intvals(:,3),'rd-','MarkerSize',15);
            errorbar(u2_intvals(:,1),u2_intvals(:,2),u2_intvals(:,3),'bd-','MarkerSize',15);
            errorbar(u2_intvals(:,1),u2_intvals(:,2)+u1_intvals(:,2),mean([u1_intvals(:,3) u2_intvals(:,3)],2),'kd-','MarkerSize',15);
            xlabel('xvalues (scale)');
            ylabel('I (um^2/s)');
            axis([-0.05 1.05 -300 200])
            xticks=([0 0.25 0.5 0.75 1]);
            yticks([-300 -300 -200 -100 0 100 200 300]);
            hold off
            
            obj.integralvalues=[u1_intvals ; u2_intvals];
        end
        function obj=calc_integral_values2(obj,listofxvals)
            %calculates and plots the integral under the u1 and u2 curves
            %of the ui vs r curves as a function of the control parameter
            %The interpretation of the values is the dissipated energy
            
            %for i=1:condip.number_of_droplets
            %        
            %    end
            
            u1_intvals=NaN(length(listofxvals),3);
            u2_intvals=NaN(length(listofxvals),3);
            u1_intvals(:,1)=listofxvals';
            u2_intvals(:,1)=listofxvals';
            
            figure
            hold on
            box on
            for p=1:obj.numberofconditions
                condip=obj.listofconditions{p};
                u1_p=condip.u1_average;
                u2_p=condip.u2_average;
                
                
                u1_int=nansum(u1_p(:,2));
                u2_int=nansum(u2_p(:,2));
                
                %error calculation from error propagation of the stdev of
                %the integral calculation
                u1_error=sqrt(nansum(u1_p(:,3).^2));
                u2_error=sqrt(nansum(u2_p(:,3).^2));
                
                %error calculation from stdev of three experiments
                u1_error=stdev(u1_p(:,2));
                u2_error=stdev(u2_p(:,2));
                
                %error calculation as lowest and highest measured values
                %and average
                                
                %u1_error=nansum(u1_p(:,3));
                %u2_error=nansum(u2_p(:,3));
                u1_intvals(p,2:3)=[u1_int u1_error];
                u2_intvals(p,2:3)=[u2_int u2_error];
            end
            %errorbar(u1_intvals(:,1),u1_intvals(:,2),u1_intvals(:,3),'rd-','MarkerSize',15);
            %errorbar(u2_intvals(:,1),u2_intvals(:,2),u2_intvals(:,3),'bd-','MarkerSize',15);
            errorbar(u2_intvals(:,1),u2_intvals(:,2)+u1_intvals(:,2),mean([u1_intvals(:,3) u2_intvals(:,3)],2),'kd-','MarkerSize',15);
            errorbar(u2_intvals(:,1),u1_intvals(:,2)-u2_intvals(:,2),mean([u1_intvals(:,3) u2_intvals(:,3)],2),'bd-','MarkerSize',15);
            xlabel('xvalues (scale)');
            ylabel('I (um^2/s)');
            axis([-0.05 1.05 -200 500])
            yticks([-300 -200 -100 0 100 200 300 400 500 600 700]);
            xticks=([0 0.25 0.5 0.75 1])
            hold off
            
            obj.integralvalues=[u1_intvals ; u2_intvals];
        end
        function obj=calc_integral_values_error(obj,listofxvals)
            %calculates and plots the integral under the u1 and u2 curves
            %of the ui vs r curves as a function of the control parameter
            %The interpretation of the values is the dissipated energy
            
            Cdrag=5.89*10^(-19); %(Ns/um)
            kT=4.11*10^-21;
            %framerate correction factor
            %!!!FOR NORMAL USE, REMOVE THIS FACTOR!!!
            fcf=1/1.2;
          
            uplus_intvals=NaN(length(listofxvals),3);
            umin_intvals=NaN(length(listofxvals),3);
            uplus_intvals(:,1)=listofxvals';
            umin_intvals(:,1)=listofxvals';
            
            figure
            hold on
            box on
            for p=1:obj.numberofconditions
                condip=obj.listofconditions{p};
                
                uplus_intvals_p=NaN(length(condip.number_of_droplets),1);
                umin_intvals_p=NaN(length(condip.number_of_droplets),1);
                                
                for i=1:condip.number_of_droplets
                    droplet_i=condip.list_of_droplets{1,i};
                    u1_pi=droplet_i.u1;
                    u2_pi=droplet_i.u2;
                    u1_int_p=nansum(u1_pi(:,2));
                    u2_int_p=nansum(u2_pi(:,2));
                    
                    uplus_intvals_p(i,1)=u1_int_p+u2_int_p;
                    umin_intvals_p(i,1)=u1_int_p-u2_int_p;
                    
                end
            
                
                uplus_int=nanmean(uplus_intvals_p(:,1));
                umin_int=nanmean(umin_intvals_p(:,1));
                
                %error calculation from stdev of 3 experiments
                uplus_intmin=min(uplus_intvals_p(:,1));
                uplus_intmax=max(uplus_intvals_p(:,1));
                umin_intmin=min(umin_intvals_p(:,1));
                umin_intmax=max(umin_intvals_p(:,1));
                
                uplus_intvals(p,2:4)=fcf*10^-4/kT*Cdrag*[uplus_int uplus_intmax uplus_intmin];
                umin_intvals(p,2:4)=fcf*10^-4/kT*Cdrag*[umin_int umin_intmax umin_intmin];
            end
            
%             for p=1:obj.numberofconditions
%                 condip=obj.listofconditions{p};
%                 u1_p=condip.u1_average;
%                 u2_p=condip.u2_average;
%                 
%                 u1_int=nansum(u1_p(:,2));
%                 u2_int=nansum(u2_p(:,2));
%                              
%                 %u1_error=nansum(u1_p(:,3));
%                 %u2_error=nansum(u2_p(:,3));
%                 u1_intvals(p,2)=[u1_int];
%                 u2_intvals(p,2)=[u2_int];
%             end
            
            %errorbar(u1_intvals(:,1),u1_intvals(:,2),u1_intvals(:,3),'rd-','MarkerSize',15);
            %errorbar(u2_intvals(:,1),u2_intvals(:,2),u2_intvals(:,3),'bd-','MarkerSize',15);
            %errorbar(u2_intvals(:,1),u2_intvals(:,2)+u1_intvals(:,2),mean([u1_intvals(:,3) u2_intvals(:,3)],2),'kd-','MarkerSize',15);
            %errorbar(u2_intvals(:,1),u1_intvals(:,2)-u2_intvals(:,2),mean([u1_intvals(:,3) u2_intvals(:,3)],2),'bd-','MarkerSize',15);
            plot(uplus_intvals(:,1),uplus_intvals(:,2),'kd-','Markersize',15);
            plot(uplus_intvals(:,1),uplus_intvals(:,3),'k.','Markersize',10);
            plot(uplus_intvals(:,1),uplus_intvals(:,4),'k.','Markersize',10);
            plot(umin_intvals(:,1),umin_intvals(:,2),'bd-','Markersize',15);
            plot(umin_intvals(:,1),umin_intvals(:,3),'b.','Markersize',10);
            plot(umin_intvals(:,1),umin_intvals(:,4),'b.','Markersize',10);
            xlabel('xvalues (scale)');
            ylabel('I (um^2/s)');
            %axis([0 0.6 -200 500])
            %yticks([-300 -200 -100 0 100 200 300 400 500 600 700]);
            %xticks=([- 0.1 0.2 0.3 0.4 0.5 0.6]);
            hold off
            
            obj.integralvalues=[uplus_intvals ; umin_intvals];
        end
        function obj=calc_integral_values_error2(obj,listofxvals)
            %calculates and plots the integral under the u1 and u2 curves
            %of the ui vs r curves as a function of the control parameter
            %The interpretation of the values is the dissipated energy
            
            Cdrag=5.89*10^(-19); %(Ns/um)
            kT=4.11*10^-21;
            %framerate correction factor
            %!!!FOR NORMAL USE, REMOVE THIS FACTOR!!!
            fcf=1/1.2;
          
            uplus_intvals=NaN(length(listofxvals),3);
            umin_intvals=NaN(length(listofxvals),3);
            uplus_intvals(:,1)=listofxvals';
            umin_intvals(:,1)=listofxvals';
            
            figure
            hold on
            box on
            for p=1:obj.numberofconditions
                condip=obj.listofconditions{p};
                
                uplus_intvals_p=NaN(length(condip.number_of_droplets),1);
                umin_intvals_p=NaN(length(condip.number_of_droplets),1);
                                
                for i=1:condip.number_of_droplets
                    droplet_i=condip.list_of_droplets{1,i};
                    u1_pi=droplet_i.u1_average;
                    u2_pi=droplet_i.u2_average;
                    u1_int_p=nansum(u1_pi(:,2));
                    u2_int_p=nansum(u2_pi(:,2));
                    
                    uplus_intvals_p(i,1)=u1_int_p+u2_int_p;
                    umin_intvals_p(i,1)=u1_int_p-u2_int_p;
                    
                end
            
                
                uplus_int=nanmean(uplus_intvals_p(:,1));
                umin_int=nanmean(umin_intvals_p(:,1));
                
                %error calculation from stdev of 3 experiments
                uplus_intmin=min(uplus_intvals_p(:,1));
                uplus_intmax=max(uplus_intvals_p(:,1));
                umin_intmin=min(umin_intvals_p(:,1));
                umin_intmax=max(umin_intvals_p(:,1));
                
                uplus_intvals(p,2:4)=2.18/1.32*fcf*10^-4/kT*Cdrag*[uplus_int uplus_intmax uplus_intmin];
                umin_intvals(p,2:4)=2.18/1.32*fcf*10^-4/kT*Cdrag*[umin_int umin_intmax umin_intmin];
            end
            
            plot(uplus_intvals(:,1),uplus_intvals(:,2),'kd-','Markersize',15);
            plot(uplus_intvals(:,1),uplus_intvals(:,3),'k.','Markersize',10);
            plot(uplus_intvals(:,1),uplus_intvals(:,4),'k.','Markersize',10);
            plot(umin_intvals(:,1),umin_intvals(:,2),'bd-','Markersize',15);
            plot(umin_intvals(:,1),umin_intvals(:,3),'b.','Markersize',10);
            plot(umin_intvals(:,1),umin_intvals(:,4),'b.','Markersize',10);
            xlabel('xvalues (scale)');
            ylabel('I (um^2/s)');
            axis([-0.05 1.05 -4 6])
            yticks([-4 -2 0 2 4 6]);
            xticks=([0 0.25 0.5 0.75 1.0]);
            hold off
            
            obj.integralvalues=[uplus_intvals ; umin_intvals];
        end
        function obj=calc_integral_values2E(obj,listofxvals)
            %calculates and plots the integral under the u1 and u2 curves
            %of the ui vs r curves as a function of the control parameter
            %The interpretation of the values is the dissipated energy
            
            Cdrag=5.89*10^(-19); %(Ns/um)
            kT=4.11*10^-21;
            %framerate correction factor
            %!!!FOR NORMAL USE, REMOVE THIS FACTOR!!!
            fcf=1/1.2;
            
            disp('Values displayed in 10^-4 kT')
            
            u1_intvals=NaN(length(listofxvals),3);
            u2_intvals=NaN(length(listofxvals),3);
            u1_intvals(:,1)=listofxvals';
            u2_intvals(:,1)=listofxvals';
            
            figure
            hold on
            box on
            for p=1:obj.numberofconditions
                condip=obj.listofconditions{p};
                u1_p=condip.u1_average;
                u2_p=condip.u2_average;
                
                
                u1_int=nansum(u1_p(:,2));
                u2_int=nansum(u2_p(:,2));
                
                u1_error=sqrt(nansum(u1_p(:,3).^2));
                u2_error=sqrt(nansum(u2_p(:,3).^2));
                %u1_error=nansum(u1_p(:,3));
                %u2_error=nansum(u2_p(:,3));
                u1_intvals(p,2:3)=[u1_int u1_error];
                u2_intvals(p,2:3)=[u2_int u2_error];
            end
            %errorbar(u1_intvals(:,1),u1_intvals(:,2),u1_intvals(:,3),'rd-','MarkerSize',15);
            %errorbar(u2_intvals(:,1),u2_intvals(:,2),u2_intvals(:,3),'bd-','MarkerSize',15);
            errorbar(u2_intvals(:,1),fcf*10^-4/kT*Cdrag*(u2_intvals(:,2)+u1_intvals(:,2)),fcf*10^-4/kT*Cdrag*mean([u1_intvals(:,3) u2_intvals(:,3)],2),'kd-','MarkerSize',15);
            errorbar(u2_intvals(:,1),fcf*10^-4/kT*Cdrag*(u1_intvals(:,2)-u2_intvals(:,2)),fcf*10^-4/kT*Cdrag*mean([u1_intvals(:,3) u2_intvals(:,3)],2),'bd-','MarkerSize',15);
            xlabel('xvalues (scale)');
            ylabel('I (um^2/s)');
            axis([-5 105 -4 8])
            xticks([0 25 50 75 100]);
            yticks([-4 -2 0 2 4 6 8]);
            hold off
            
            obj.integralvalues=[u1_intvals ; u2_intvals];
        end
        function plot_fit_lin(obj)
            %plot the approach and escape speed of droplets on a lin lin
            %scale U vs r with fits of shape U(r) = A r^-2 + B
            cmap=winter(obj.numberofconditions);
            
            %plot figure 1
            figure
            hold on
            box on
            title('u1')
            for p=1:obj.numberofconditions
                
                condip=obj.listofconditions{p};
                
                %obtain input matrices
                u1_av=condip.u1_average;
                                               
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
                u1_clip=u1_av(find(u1_av(:,1)==rmin):find(u1_av(:,1)==rmax),:);
                u1_clip(isnan(u1_clip(:,2)),:)=[];
                
                %fit data
                [C12,S12]=polyfit(u1_clip(:,1).^(-2),u1_clip(:,2),1);
                                
                %plot data u1
                x=u1_clip(:,1);
                y=u1_clip(:,2);
                error=u1_clip(:,3);
                x1=linspace(rmin-3,rmax+3,100);
                y1=C12(1)./x1.^2;
                y1=y1(:)-y1(100);
                errorbar(x,y,error,'sq','Color',cmap(p,:))
                plot(x1,y1,'-','Color',cmap(p,:));
                
                
            end
            xlabel('r (um)');
            ylabel('v (um/s)');
            axis([0 180 -20 20])
            yticks([-20 -10 0 10 20])
            %xticks([80 100 120 140 160 180])
            xticks=([0 50 100 150]);
            hold off
            
            %plot figure 2
            figure
            hold on
            box on
            title('u2')
            for p=1:obj.numberofconditions
                
                condip=obj.listofconditions{p};
                
                %obtain input matrices
                u2_av=condip.u2_average;
                u1_av=condip.u1_average;
                
                %make sure enough arguments are provided
                if nargin<3
                    upperlim=max(u2_av(:,1));
                end
                if nargin<2
                    lowerlim=min(u2_av(:,1));
                end
                
                %find right range so that the range falls within both the
                %user provided limits and the measured distance range
                rmin=max([lowerlim min(u2_av(:,1))]);
                rmax=min([upperlim max(u2_av(:,1))]);
                
                %clip data to right range
                u2_clip=u2_av(find(u2_av(:,1)==rmin):find(u2_av(:,1)==rmax),:);
                u2_clip(isnan(u2_clip(:,2)),:)=[];
                
                %fit data
                [C21,S21]=polyfit(u2_clip(:,1).^(-2),u2_clip(:,2),1);
                             
                %plot data u2
                x=u2_clip(:,1);
                y=u2_clip(:,2);
                error=u2_clip(:,3);
                x1=linspace(rmin-3,rmax+3,100);
                y1=C21(1)./x1.^2;
                y1=y1(:)-y1(100);
                errorbar(x,y,error,'sq','Color',cmap(p,:))
                plot(x1,y1,'-','Color',cmap(p,:));
                
                
            end
            xlabel('r^ (um)');
            ylabel('v (um/s)');
            axis([0 180 -20 20])
            yticks([-20 -10 0 10 20])
            %xticks([80 100 120 140 160 180])
            xticks=([0 50 100 150]);
            hold off
        end
        function plot_ui_all_scaled(obj,listofscales)
            %plots the speed of the chaser and chasee as function of
            %interparticle distance for all conditions
            %divides y value by scaling factor indicated in listofscales
            %listofscales needs to be of form [sc1 sc2 sc3...] with length
            %equal to the number of conditions
            
            %get colormap
            cmap=winter(obj.numberofconditions);
            
            %make figure u1
            figure
            hold on
            box on
            title('U1')
            for p=1:obj.numberofconditions
                %Get values to plot
                condip=obj.listofconditions{p};
                u1_p=condip.u1_average;
                x=u1_p(:,1);
                y=u1_p(:,2)./listofscales(p);
                error=u1_p(:,3)./listofscales(p);
                %plot values
                errorbar(x,y,error,'sq','Color',cmap(p,:))
            end
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            hold off
            
            %make figure u2
            figure
            hold on
            box on
            title('U2')
            for p=1:obj.numberofconditions
                %Get values to plot
                condip=obj.listofconditions{p};
                u2_p=condip.u2_average;
                x=u2_p(:,1);
                y=u2_p(:,2)./listofscales(p);
                error=u2_p(:,3)./listofscales(p);
                %plot values
                errorbar(x,y,error,'sq','Color',cmap(p,:))
            end
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            hold off
        end
        function plot_ui_all_scaled_loglog(obj,listofscales)
            %plots the speed of the chaser and chasee as function of
            %interparticle distance for all conditions
            %divides y value by scaling factor indicated in listofscales
            %listofscales needs to be of form [sc1 sc2 sc3...] with length
            %equal to the number of conditions
            
            %get colormap
            cmap=winter(obj.numberofconditions);
            
            %make figure u1
            figure
            hold on
            box on
            title('U1')
            for p=1:obj.numberofconditions
                %Get values to plot
                condip=obj.listofconditions{p};
                u1_p=condip.u1_average;
                x=u1_p(:,1);
                y=abs(u1_p(:,2)./listofscales(p));
                error=u1_p(:,3)./listofscales(p);
                %plot values
                errorbar(x,y,error,'sq','Color',cmap(p,:))
                y(1)
                y(end)
                plot(x,0.6*10^4*x.^-2,'r-')
                plot(x,10^21*x.^-10,'b-')
            end
            set(gca, 'YScale', 'log')
            set(gca, 'XScale', 'log')
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            hold off
            
            %make figure u2
            figure
            hold on
            box on
            title('U2')
            for p=1:obj.numberofconditions
                %Get values to plot
                condip=obj.listofconditions{p};
                u2_p=condip.u2_average;
                x=u2_p(:,1);
                y=abs(u2_p(:,2)./listofscales(p));
                error=u2_p(:,3)./listofscales(p);
                %plot values
                errorbar(x,y,error,'sq','Color',cmap(p,:))
                plot(x,0.6*10^4*x.^-2,'r-')
                plot(x,10^21*x.^-10,'b-')
            end
            set(gca, 'YScale', 'log')
            set(gca, 'XScale', 'log')
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            hold off
        end
        function plot_u_rel_all(obj)
            %plots the speed of the chaser and chasee as function of
            %interparticle distance for all conditions
            
            %get colormap
            cmap=winter(obj.numberofconditions);
            
            %make figure u_rel
            figure
            hold on
            box on
            title('u_rel')
            for p=1:obj.numberofconditions
                %Get values to plot
                condip=obj.listofconditions{p}
                p
                u1_p=condip.u_rel_average;
                x=u1_p(:,1);
                y=u1_p(:,2);
                error=u1_p(:,3);
                %plot values
                errorbar(x,y,error,'sq','Color',cmap(p,:))
            end
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            hold off
            
        end
        function plot_u_sum_all(obj)
            %plots the speed of the chaser and chasee as function of
            %interparticle distance for all conditions
            
            %get colormap
            cmap=winter(obj.numberofconditions);
            
            %make figure u sum
            figure
            hold on
            box on
            title('u sum')
            for p=1:obj.numberofconditions
                %Get values to plot
                condip=obj.listofconditions{p};
                u1_p=condip.u_sum_average;
                x=u1_p(:,1);
                y=u1_p(:,2);
                error=u1_p(:,3);
                %plot values
                errorbar(x,y,error,'sq','Color',cmap(p,:))
            end
            xlabel('r (\mum)');
            ylabel('v (\mum/s)');
            hold off
            
        end
        function plot_ui_rs_all(obj)
            %plots the speed of the chaser and chasee as function of
            %interparticle distance to the power -2 for all conditions
            %get colormap
            cmap=winter(obj.numberofconditions);
            
            %make figure u1
            figure
            hold on
            box on
            title('U1')
            for p=1:obj.numberofconditions
                %Get values to plot
                condip=obj.listofconditions{p};
                u1_p=condip.u1_average;
                x=u1_p(:,1).^(-2);
                y=u1_p(:,2);
                error=u1_p(:,3);
                %plot values
                errorbar(x,y,error,'sq','Color',cmap(p,:))
            end
            xlabel('r^{-2} (\mum^{-2})');
            ylabel('v (\mum/s)');
            hold off
            
            %make figure u2
            figure
            hold on
            box on
            title('U2')
            for p=1:obj.numberofconditions
                %Get values to plot
                condip=obj.listofconditions{p};
                u2_p=condip.u2_average;
                x=u2_p(:,1).^(-2);
                y=u2_p(:,2);
                error=u2_p(:,3);
                %plot values
                errorbar(x,y,error,'sq','Color',cmap(p,:))
            end
            xlabel('r^{-2} (\mum^{-2})');
            ylabel('v (\mum/s)');
            hold off
            
        end
        function plot_fits(obj,lowerlim,upperlim)
            %plots the fits that give the products C12 and C21 for all conditions
            
            cmap=winter(obj.numberofconditions);
            
            %plot figure 1
            figure
            hold on
            box on
            title('u1')
            for p=1:obj.numberofconditions
                
                condip=obj.listofconditions{p};
                
                %obtain input matrices
                u1_av=condip.u1_average;
                                               
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
                u1_clip=u1_av(find(u1_av(:,1)==rmin):find(u1_av(:,1)==rmax),:);
                u1_clip(isnan(u1_clip(:,2)),:)=[];
                
                %fit data
                [C12,S12]=polyfit(u1_clip(:,1).^(-2),u1_clip(:,2),1);
                                
                %plot data u1
                x=u1_clip(:,1).^(-2);
                y=u1_clip(:,2);
                error=u1_clip(:,3);
                x1=linspace(rmin-3,rmax+3,100).^(-2);
                y1=polyval(C12,x1);
                errorbar(x,y,error,'sq','Color',cmap(p,:))
                plot(x1,y1,'-','Color',cmap(p,:));
                xlabel('r^{-2} (\mum^{-2})');
                ylabel('v (\mum/s)');
                
            end
            xlabel('r^{-2} (\mum^{-2})');
            ylabel('v (\mum/s)');
            axis([0.2*10^-4 1.2*10^-4 -20 20])
            yticks([-20 -10 0 10 20])
            xticks(10^-5*[2 4 6 8 10 12])
            hold off
            
            %plot figure 2
            figure
            hold on
            box on
            title('u2')
            for p=1:obj.numberofconditions
                
                condip=obj.listofconditions{p};
                
                %obtain input matrices
                u2_av=condip.u2_average;
                u1_av=condip.u1_average;
                
                %make sure enough arguments are provided
                if nargin<3
                    upperlim=max(u2_av(:,1));
                end
                if nargin<2
                    lowerlim=min(u2_av(:,1));
                end
                
                %find right range so that the range falls within both the
                %user provided limits and the measured distance range
                rmin=max([lowerlim min(u2_av(:,1))]);
                rmax=min([upperlim max(u2_av(:,1))]);
                
                %clip data to right range
                u2_clip=u2_av(find(u2_av(:,1)==rmin):find(u2_av(:,1)==rmax),:);
                u2_clip(isnan(u2_clip(:,2)),:)=[];
                
                %fit data
                [C21,S21]=polyfit(u2_clip(:,1).^(-2),u2_clip(:,2),1);
                             
                %plot data u2
                x=u2_clip(:,1).^(-2);
                y=u2_clip(:,2);
                error=u2_clip(:,3);
                x1=linspace(rmin-3,rmax+3,100).^(-2);
                y1=polyval(C21,x1);
                errorbar(x,y,error,'sq','Color',cmap(p,:))
                plot(x1,y1,'-','Color',cmap(p,:));
                xlabel('r^{-2} (\mum^{-2})');
                ylabel('v (\mum/s)');
                
            end
            xlabel('r^{-2} (\mum^{-2})');
            ylabel('v (\mum/s)');
            axis([0.2*10^-4 1.2*10^-4 -20 20])
            yticks([-20 -10 0 10 20])
            xticks(10^-5*[2 4 6 8 10 12])
            hold off
        end
            
    end %end methods
end %end class