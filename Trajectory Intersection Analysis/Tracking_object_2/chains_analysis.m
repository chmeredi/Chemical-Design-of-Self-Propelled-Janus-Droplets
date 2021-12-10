classdef chains_analysis
    %Class designed to average properties of object of type single chain
    %All chains have to be of the same length for this code to work (or
    %make sense)
    
    properties
        scale
        framerate
        chainlength
        listofchains
        numberofchains
        tensor_average
        corrs
        tensor_stdev
        endvector_all
        contourlength_all
        D_trans
        R_end
        NOF %largest number of frames in list of objects
    end
    methods
        function obj=chains_analysis(listofsinglechains)
            %list of single chains must be a 1xn cell {obj1 obj2...}
            %of objects of type single chain
            
            
            obj.listofchains={};
            obj1=listofsinglechains{1};
            NOP_1=obj1.NOP;
            
            for p=1:length(listofsinglechains)
                chain_p=listofsinglechains{p};
                obj.listofchains{p}=chain_p;
                NOP_temp=chain_p.NOP;
                if NOP_temp~=NOP_1
                    error('the chain length is not the same for all the objects you gave as input')
                end
            end
            obj.numberofchains=length(listofsinglechains);
            obj.scale=obj1.scale;
            obj.framerate=obj1.framerate;
            obj.chainlength=obj1.chainlength;
            
            maxval=1;
            for p=1:length(obj.listofchains)
                obj_p=obj.listofchains{p};
                NOF_p=obj_p.NOF;
                if NOF_p>maxval
                    maxval=NOF_p;
                end
            end
            obj.NOF=maxval;
        end
        function obj=add_chain(obj,chain)
            %adds an instance of dropletpair to the object
            
            listtemp=obj.listofchains;
            listtemp{end+1}=chain;
            obj.listofchains=listtemp;
            obj.numberofchains=obj.numberofchains+1;
            
        end
        function obj=remove_chain(obj,p)
            %removes particle with index pfrom the list of dropletpairs
            
            if p<=obj.numberofchains
                listnew={};
                count=0;
                listtemp=obj.listofchains;
                for q=1:obj.numberofchains
                    if q~=p
                        count=count+1;
                        listnew{count}=listtemp{q};
                    end
                end
            else
                disp(['There is no chain with index ' num2str(p) '.'])
            end
            obj.listofchains=listnew;
            obj.numberofchains=obj.numberofchains-1;
        end
        function obj=average_diffusion_tensor(obj)
            %averages the diffusion tensor of all enclosed objects
            allNOF=0;
            obj1=obj.listofchains{1};
            tensorsize=size(obj1.diffs);
            tensor_all=NaN(tensorsize(1),tensorsize(2),obj.numberofchains);
            tensor_allstd=NaN(tensorsize(1),tensorsize(2),obj.numberofchains);
            w=NaN(tensorsize(1),tensorsize(2),obj.numberofchains);
            for p=1:obj.numberofchains
                obj_p=obj.listofchains{p};
                if isempty(obj_p.diffs)
                    obj_p=calc_diff_tensor(obj_p);
                end
                tensor_p=obj_p.diffs;
                tensor_p=cell2mat(tensor_p)*obj_p.NOF;
                tensor_all(:,:,p)=tensor_p;
                tensor_allstd(:,:,p)=cell2mat(obj_p.diffs);
                w(:,:,p)=obj_p.NOF;
                allNOF=allNOF+obj_p.NOF;
            end
            for i=1:length(tensor_p(:,1))
                for j=1:length(tensor_p(1,:))
                    obj.tensor_stdev(i,j)=sqrt(nanvar(permute(tensor_allstd(i,j,:),[3 2 1]),permute(w(i,j,:),[3 2 1])));
                end
            end
            obj.tensor_average=nansum(tensor_all,3);
            obj.tensor_average=obj.tensor_average./allNOF;
        end
        function obj=average_corrs(obj)
            %Average all correlations functions in the diffusion tensor
            
            obj1=obj.listofchains{1};
            ni=length(obj1.diffs);
            corrs_av=cell(ni,ni);
            
            maxval=floor(obj.NOF/10);
                       
            for i=1:ni
                for j=1:ni
                    Cij_all=NaN(maxval,1);
                    for p=1:obj.numberofchains
                        obj_p=obj.listofchains{p};
                        corr_p=obj_p.corrs;
                        Cij=corr_p{i,j};
                        Cij_all(1:length(Cij(:,2)),p)=Cij(:,2);
                    end
                    Cij_av=nanmean(Cij_all,2);
                    Cij_std=nanstd(Cij_all');
                    Cij_std=Cij_std';
                    corrs_av{i,j}=[linspace(0,length(Cij_av(:,1))-1,length(Cij_av(:,1)))' Cij_av Cij_std];
                end
            end
            
            obj.corrs=corrs_av;
            
        end
        function plot_corr_funcs(obj)
            %plots the x and y MSD graph in one graph
            %plots the rotational diffusion graph
            %plots the internal bond angle graph
            
            if isempty(obj.corrs)
                obj=average_corrs(obj);
            end
            if isempty(obj.tensor_average)
                obj=average_diffusion_tensor(obj);
            end
            
            %x and y MSD graph
            Cx=obj.corrs{1,1};
            Cy=obj.corrs{2,2};
            Dx=obj.tensor_average(1,1)
            Dy=obj.tensor_average(2,2)
            figure
            hold on
            box on
            errorbar(Cx(:,1)/obj.framerate,Cx(:,2)*obj.scale^2,Cx(:,3)*obj.scale^2,'bo');
            errorbar(Cy(:,1)/obj.framerate,Cy(:,2)*obj.scale^2,Cy(:,3)*obj.scale^2,'ro');
            plot(Cx(:,1)/obj.framerate,Cx(:,1)/obj.framerate*Dx,'b-');
            plot(Cy(:,1)/obj.framerate,Cy(:,1)/obj.framerate*Dy,'r-');
            xlabel('dt (s)');
            ylabel('MSD (\um^2)');
            hold off
            
            %rot diff graph
            Cr=obj.corrs{3,3};
            Dr=obj.tensor_average(3,3);
            figure
            hold on
            box on
            errorbar(Cr(:,1)/obj.framerate,Cr(:,2),Cr(:,3),'ko');
            plot(Cr(:,1)/obj.framerate,Cr(:,1)/obj.framerate*Dr,'r-');
            xlabel('dt (s)');
            ylabel('MSD (rad^2)');
            hold off
            
            %alpha's
            for p=4:obj.chainlength+1
                Ca=obj.corrs{p,p};
                Da=obj.tensor_average(p,p);
                figure
                hold on
                box on
                title(['Bond angle ' num2str(p-3) '.']);
                errorbar(Ca(:,1)/obj.framerate,Ca(:,2),Ca(:,3),'ko');
                plot(Ca(:,1)/obj.framerate,Ca(:,1)/obj.framerate*Da,'r-');
                xlabel('dt (s)');
                ylabel('MSD (rad^2)');
                hold off
            end
            
            
        end
        function plot_all_correlations(obj)
        %Plot the average correlation functions
            NOP_1=obj.chainlength;
            average_corrs=obj.corrs;
            f = figure(10);
            p = uipanel('Parent',f,'BorderType','none','BackgroundColor','white');
            p.TitlePosition = 'centertop';
            p.FontSize = 12;
            p.FontWeight = 'bold';
            for i=1:NOP_1+1
                for j=1:NOP_1+1
                    Cij=average_corrs{i,j};
                    hold on
                    box on
                    axi=subplot(NOP_1+1,NOP_1+1,(j-1)*(NOP_1+1)+i,'Parent',p);
                    title(['C' num2str(i) num2str(j)])
                    
                    if (i==1 || i==2) && (j==1 || j==2)
                        errorbar(axi,Cij(:,1)/obj.framerate,Cij(:,2)*obj.scale^2,Cij(:,3)*obj.scale^2,'ko');
                        ylabel('C(t) (\mum^2)');
                    elseif i==1 || i==2 || j==1 || j==2
                        errorbar(axi,Cij(:,1)/obj.framerate,Cij(:,2)*obj.scale,Cij(:,3)*obj.scale,'ko');
                        ylabel('C(t) (\mum.rad)');
                    else
                        errorbar(axi,Cij(:,1)/obj.framerate,Cij(:,2),Cij(:,3),'ko');
                        ylabel('C(t) (rad^2)');
                    end
                    xlabel('dt (s)')
                    if i==j
                        axis([0 10 -1 1]);
                    else
                        axis([0 10 -1 1]);
                    end
                    hold off
                    
                end
            end
        end
        function obj=plot_MSD_all(obj,far)
            %plots the mean squared displacement for all enclosed objects
            
            if nargin<2
                far=500;
            end
            
            slope=NaN(obj.numberofchains,1);
            slopestd=NaN(obj.numberofchains,1);
            slopelog=NaN(obj.numberofchains,1);
            slopelogstd=NaN(obj.numberofchains,1);
            w=NaN(obj.numberofchains,1);
            NOFsum=0;
            
            fps=obj.framerate;
            scl=obj.scale;
            
            maxval=1;
            for p=1:length(obj.listofchains)
                obj_p=obj.listofchains{p};
                NOF_p=obj_p.NOF;
                if NOF_p>maxval
                    maxval=NOF_p;
                end
            end
            maxval=ceil(maxval/10);
            
            figure
            hold on
            box on
            cmap=jet(obj.numberofchains);
            for p=1:obj.numberofchains
                obj_p=obj.listofchains{p};
                obj_p=find_center_of_mass(obj_p);
                cmass_p=store_cmass(obj_p);
                cmass_p=calc_MSD(cmass_p);
                MSD_p=cmass_p.MSD{1,2};
                plot(MSD_p(:,1)*obj.scale^2,MSD_p(:,2)*obj.framerate,'.','Color',cmap(p,:));
                dt=(MSD_p(2,1)-MSD_p(1,1));
                MSD_log=log10(MSD_p);
                slope_temp=(MSD_p(2:end,2)-MSD_p(1,2))./(linspace(1,length(MSD_p(2:end,2)),length(MSD_p(2:end,2)))');
                slope(p,1)=nanmean(slope_temp)*obj_p.NOF;
                slopestd(p,1)=nanmean(slope_temp);
                slopelog_temp=(MSD_log(3:end,2)-MSD_log(2,2))./log10(linspace(2,length(MSD_log(3:end,2)),length(MSD_log(3:end,2)))');
                %slopelog(p,1)=nanmean(((MSD_log(3:end,2))-MSD_log(2:end-1,2))./(MSD_log(3:end,1)-MSD_log(2:end-1,1)))
                slopelog(p,1)=nanmean(slopelog_temp)*obj_p.NOF;
                slopelogstd(p,1)=nanmean(slopelog_temp);
                w(p,1)=obj_p.NOF;
                NOFsum=NOFsum+obj_p.NOF;
                if p==1
                    MSD_all=NaN(maxval,obj.numberofchains);  
                end
                MSD_all(1:length(MSD_p(:,2)),p)=MSD_p(:,2);
            end
            D_trans_temp=scl^2*fps*[nansum(slope(:,1))/4/NOFsum sqrt(nanvar(slopestd(:,1),w))/4];
            obj.D_trans(1:2,1)=D_trans_temp;
            slopeval(1:2,1)=[nansum(slopelog(:,1))/NOFsum sqrt(nanvar(slopelogstd(:,1),w))];
            hold off
            
            %calc average MSD
            MSD_mean(:,1)=linspace(0,maxval-1,maxval)';
            MSD_mean(:,2)=nanmean(MSD_all(:,:),2);
            std_temp=MSD_all(:,:)';
            std_temp=nanstd(std_temp(:,:),1);
            MSD_mean(:,3)=std_temp;
            
            %plot the average MSD
            figure
            hold on
            box on
            errorbar(MSD_mean(1:far,1)/fps,MSD_mean(1:far,2)*scl^2,MSD_mean(1:far,3)*scl^2,'ko');
            plot(MSD_mean(1:far,1)/fps,4*obj.D_trans(1,1)*MSD_mean(1:far,1)/fps,'r-');
            xlabel('dt (s)')
            ylabel('MSD (\um^2)')
            hold off
            
            %plot the average MSD on log scale
            figure
            hold on
            box on
            h = errorbar(MSD_mean(1:far,1)/fps,MSD_mean(1:far,2)*scl^2,MSD_mean(1:far,3)*scl^2,'ko');
            h = plot(linspace(0,far-1,far)/fps,(linspace(0,far-1,far)/fps).^slopeval(1,1)*obj.D_trans(1,1)*4,'r-');
            set(gca,'YScale','log');
            set(gca,'XScale','log');
            xlabel('dt (s)')
            ylabel('MSD (\um^2)')
            hold off
           
            disp(['Average diffusion coefficient is: ' num2str(obj.D_trans(1,1)) ' um^2/s.'])
            disp(['Standard deviation is: ' num2str(obj.D_trans(2,1)) ' um^2/s.'])
            
            disp(['Average slope of loglogplot is: ' num2str(slopeval(1,1)) '.'])
            disp(['Standard deviation is: ' num2str(slopeval(2,1)) '.'])
        end
        function obj=plot_endvector_histogram_all(obj)
            %plots the average histogram of endvectors
            
            obj.endvector_all=[];
            obj.contourlength_all=[];
            for p=1:obj.numberofchains
                obj_p=obj.listofchains{p};
                if isempty(obj_p.endvector)
                    obj_p=calc_endvector(obj_p);
                end
                if isempty(obj_p.contourlength)
                    obj_p=find_contourlength(obj_p);
                end
                endvector_p=obj_p.endvector;
                endvector_p(:,1)=endvector_p(:,1)*obj_p.scale;
                obj.endvector_all=[obj.endvector_all ; endvector_p];
                obj.contourlength_all=[obj.contourlength_all ; obj_p.contourlength*obj_p.scale];
            end
            
            figure
            hold on
            box on
            legend('show')
            histogram(obj.endvector_all(:,1),'DisplayName','End-to-end distance')
            histogram(obj.contourlength_all(:,1),'DisplayName','Contour length')
            xlabel('length (\mum)','FontSize',18)
            ylabel('Number','FontSize',18)
            hold off
            
            obj.R_end(1:2,1)=[nanmean(obj.endvector_all(:,1)) nanstd(obj.endvector_all(:,1))];
            
            disp(['Average end-to-end distance is: ' num2str(nanmean(obj.endvector_all(:,1))) ' um.'])
            disp(['Standard deviation is: ' num2str(nanstd(obj.endvector_all(:,1))) ' um.'])
            
        end
        function plot_DxDy_trend_angle_all(obj,anglerange,angleav)
            %WORKS ONLY FOR TRIMERS
            %plots Dx, Dy and Dtrans as function of time
            %interval used for the measurement
            %
            %anglerange is width of accepted angles starting from average
            %angle that is accepted in radians. 2pi/3 means all angles.
            %
            %angleav is the average angle around which this is calculated
            %0 is stretched, 2 pi/3 is a closed triangle
                        
            if nargin<2
                anglerange=0.1;
            end
            if nargin<3
                angleav=0;
            end
            
            tmax=obj.NOF;
            dtmax=floor(tmax/10);
            
            Dt_p=NaN(dtmax,obj.numberofchains);
            Dx_p=NaN(dtmax,obj.numberofchains);
            Dy_p=NaN(dtmax,obj.numberofchains);
            for p=1:obj.numberofchains
                Dxy=NaN(dtmax,3);
                obj_p=obj.listofchains{p};
                tmax_temp=obj_p.NOF;
                %for each time interval
                for dt=1:dtmax
                    %for each initial position
                    Cxy=NaN(tmax-dt,2);
                    for t=1:dt:tmax_temp-dt
                        anglest=obj_p.angles{1,2};
                        angle_i=anglest(t,1);
                        %if the initial angle is in the right range
                        if abs(angle_i)>abs(angleav)-anglerange && abs(angle_i)<abs(angleav)+anglerange
                            %the orientation angle in the intial position is
                            angle=-obj_p.endvector(t,2);
                            %gives a rotation matrix R
                            R=[cos(angle) -sin(angle) ; sin(angle) cos(angle)];
                            %get the trajectory for the object from time t onward
                            %and translate such that position at dt=0 is at origin
                            traj=obj_p.tr_mass(t+dt,1:2)-obj_p.tr_mass(t,1:2);
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

                temp1=store_cmass(obj_p);
                temp1=calc_MSD(temp1);
                MSD_t=temp1.MSD{1,2};
                Dt_1=MSD_t(2:end,1);
                Dt_p(1:length(Dt_1(:,1)),p)=(1/4)*MSD_t(2:end,2)./MSD_t(2:end,1);
                Dx_p(1:length(Dxy(:,2)),p)=Dxy(:,2);
                Dy_p(1:length(Dxy(:,3)),p)=Dxy(:,3);
            end
            Dt=[linspace(1,length(Dt_p(:,1)),length(Dt_p(:,1)))' nanmean(Dt_p(:,:),2)];
            Dx=nanmean(Dx_p(:,:),2);
            Dy=nanmean(Dy_p(:,:),2);
            
            figure
            box on
            semilogx(Dt(1:dtmax,1)/obj.framerate,Dt(1:dtmax,2)*obj.scale^2*obj.framerate,'k-',...
            Dxy(1:dtmax,1)/obj.framerate,Dy(1:dtmax,1)*obj.scale^2*obj.framerate,'ro-',...
            Dxy(1:dtmax,1)/obj.framerate,Dx(1:dtmax,1)*obj.scale^2*obj.framerate,'bo-')
            
        
           
        end
        function plot_DxDy_trend_all(obj,tmax)
            %plots Dx, Dy and Dtrans as function of time
            %interval used for the measurement
            
            if nargin<2 || tmax>obj.NOF
                tmax=obj.NOF;
            end
            
            dtmax=floor(tmax/10);
            
            Dt_p=NaN(dtmax,obj.numberofchains);
            Dx_p=NaN(dtmax,obj.numberofchains);
            Dy_p=NaN(dtmax,obj.numberofchains);
            for p=1:obj.numberofchains
                Dxy=NaN(dtmax,3);
                obj_p=obj.listofchains{p};
                tmax_temp=obj_p.NOF;
                %for each time interval
                for dt=1:dtmax
                    %for each initial position
                    Cxy=NaN(tmax-dt,2);
                    for t=1:dt:tmax_temp-dt
                        %the orientation angle in the intial positiston is
                        angle=-obj_p.endvector(t,2);
                        %gives a rotation matrix R
                        R=[cos(angle) -sin(angle) ; sin(angle) cos(angle)];
                        %get the trajectory for the object from time t onward
                        %and translate such that position at dt=0 is at origin
                        traj=obj_p.tr_mass(t+dt,1:2)-obj_p.tr_mass(t,1:2);
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

                temp1=store_cmass(obj_p);
                temp1=calc_MSD(temp1);
                MSD_t=temp1.MSD{1,2};
                Dt_1=MSD_t(2:end,1);
                Dt_p(1:length(Dt_1(:,1)),p)=(1/4)*MSD_t(2:end,2)./MSD_t(2:end,1);
                Dx_p(1:length(Dxy(:,2)),p)=Dxy(:,2);
                Dy_p(1:length(Dxy(:,3)),p)=Dxy(:,3);
            end
            Dt=[linspace(1,length(Dt_p(:,1)),length(Dt_p(:,1)))' nanmean(Dt_p(:,:),2)];
            Dx=nanmean(Dx_p(:,:),2);
            Dy=nanmean(Dy_p(:,:),2);
            
            figure
            box on
            semilogx(Dt(1:dtmax,1)/obj.framerate,Dt(1:dtmax,2)*obj.scale^2*obj.framerate,'k-',...
            Dxy(1:dtmax,1)/obj.framerate,Dy(1:dtmax,1)*obj.scale^2*obj.framerate,'ro-',...
            Dxy(1:dtmax,1)/obj.framerate,Dx(1:dtmax,1)*obj.scale^2*obj.framerate,'bo-')
        
           
        end
    end%end methods
end%end class