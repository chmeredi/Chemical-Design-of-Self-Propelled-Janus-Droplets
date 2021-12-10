classdef droplets < particles
    %Droplets class is a sublcass of particles
    %Extra properties
    %   - tr_red
    %Extra methods
    %   - find_reduced_trajectory(obj)
    %   - plot_reduced_trajectory(obj,cutoffvalue)
    %   - find_pos_hough (gebruikt hough transform ipv gauss peakfind)
    %   - find_neighbours(obj,cutoffvalue)
    %   - contour(obj,npoints,factor)
    %   - output=split_to_pair(obj)
    %   - output=split_to_pair2(obj)
    %   - obj=pick_selection(obj,n)
    %Created on 01-09-17 by Pepijn Moerman
    %Last modified: 04-09-17
    
    properties
        contourfile
        %selection1 % User picked subselection of particles present in first frame. List of ID's in obj.tr;
        %selection2 % User picked subselection of particles present in first frame. List of ID's in obj.tr;
    end
    methods
        function obj=droplets(fname,parameters,identity)
            %Droplets needs to be an object of type particles
            %Identity indicates whether or not users is requested to find
            %positions and fill in scale.
            
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
            
            obj=obj@particles(fname,parameters);
            
            if identity=='y'
                obj.scale=get_numerical_input('What is the scale of the image in micrometers per pixel ?');
                obj.framerate=get_numerical_input('What is the framerate of the image? ');
                
                message='Do you want to analyze the data? Yes (y) or no (n): ';
                acceptables=['y' 'n'];
                answer=get_textual_input(message,acceptables);
                answer1='n';
                answer2='o';
                while answer=='y'
                    if answer1=='n' || answer2=='n'
                        if answer2=='o'
                            obj=find_pos_hough(obj,initialize);
                            message='Are you happy with the particle finding? Yes (y) or no (n): ';
                            answer1=get_textual_input(message,acceptables);
                        end
                        if answer1=='y' && answer2~='y'
                            proceed=1;
                            try obj=track_obj(obj,initialize);
                            catch
                                disp('Trackingcode broke. Try different parameters or go back to particlefinding.');
                                proceed=0;
                            end
                            if proceed==1
                                if ~isempty(obj.tr)
                                    plot_tracks(obj);
                                    message='Are you happy with the particle tracking? Yes (y) or no (n): ';
                                    answer2=get_textual_input(message,acceptables);
                                else
                                    answer='n';
                                    disp('Tracking failed. Try different parameters');
                                end
                            end
                            
                        end
                        if answer2=='y'
                            answer='n';
                        end
                        disp('Basic image analysis is done. Use more specific functions of the object to proceed.')
                        disp('A list of functions and properties of the object is indicated at the top of the code.')
                    end
                end
                
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
        function output=split_to_pair(obj,list)
            %Ask user to point out two particles in the first frame of a movie and outputs the trace of
            %particles as an object of type dropletpair
            
            %Asserts that the trace exits
            if isempty(obj.tr)
                obj.tr=track(obj);
            end
            
            %Asserts that the image exists
            if isempty(obj.im)
                error('Object contains no image. This way you cannot pick the right particles.');
            end
            
            trace=obj.tr;
            if obj.NOP>2 && nargin<2
                %Asks user to pick two particles
                disp('Point out the two particles you want to split to a pair: ')
                picture=obj.im;
                %picture=rgb2gray(picture);
                figure
                colormap('gray'),imagesc(picture)
                [x,y]=ginput(2);
                
                %Calculate distances of all particles in first frame to user
                %input points
                distances=NaN(obj.NOP,2);
                tr1=trace(find(trace(:,3)==1),:);
                for p=1:obj.NOP
                    trp=tr1(find(tr1(:,4)==p),:);
                    distances(p,1)=sqrt((trp(1,1)-x(1))^2+(trp(1,2)-y(1))^2);
                    distances(p,2)=sqrt((trp(1,1)-x(2))^2+(trp(1,2)-y(2))^2);
                end
                %Find minimal distances and corresponding particle IDs
                min1=nanmin(distances(:,1));
                min2=nanmin(distances(:,2));
                ID1=find(distances(:,1)==min1);
                ID2=find(distances(:,2)==min2);
                
                %Get trace of particle IDS
                tr_1=trace(find(trace(:,4)==ID1),:);
                tr_2=trace(find(trace(:,4)==ID2),:);
                tracenew=[tr_1(:,1:3) ones(length(tr_1(:,1))) ; tr_2(:,1:3) 2*ones(length(tr_1(:,1)))];
                posnew=[tr_1(:,1:3) ; tr_2(:,1:3)];
            elseif obj.NOP==2
                tracenew=trace;
                posnew=trace(:,1:3);
            elseif nargin==2
                if length(list)>2
                    list=list(1:2);
                end
                ID1=list(1);
                ID2=list(2);
                %Get trace of particle IDS
                tr_1=trace(find(trace(:,4)==ID1),:);
                tr_2=trace(find(trace(:,4)==ID2),:);
                tracenew=[tr_1(:,1:3) ones(length(tr_1(:,1))) ; tr_2(:,1:3) 2*ones(length(tr_1(:,1)))];
                posnew=[tr_1(:,1:3) ; tr_2(:,1:3)];
            else
                error('You did not give enough input arguments or your object contains fewer than 2 droplets');
            end
            
            %Create dropletpair object as output, with the right traces and
            %the parameters from this object
            newname=[obj.foldername obj.filename '.tif'];
            params=get_parameters(obj);
            output=dropletpair(newname,params,'n');
            output.scale=obj.scale;
            output.framerate=obj.framerate;
            output=add_pos(output,posnew);
            output=add_trace(output,tracenew);
            output.im=obj.im;
            output.NOP=2;
        end
        function output=split_to_pair2(obj)
            %Ask user to point out two particles in the first frame of a movie and outputs the trace of
            %particles as an object of type dropletpair2
            %IMPORTANT: CLICKING ORDER MATTERS!
            %ALWAYS FIRST SELECT CHASEE, THEN CHASER!
            
            %Asserts that the trace exits
            if isempty(obj.tr)
                obj.tr=track(obj);
            end
            
            %Asserts that the image exists
            if isempty(obj.im)
                error('Object contains no image. This way you cannot pick the right particles.');
            end
            
            trace=obj.tr;
            if obj.NOP>1
                %Asks user to pick two particles
                disp('Point out the two particles you want to split to a pair: ')
                picture=obj.im;
                %picture=rgb2gray(picture);
                figure
                colormap('gray'),imagesc(picture)
                [x y]=ginput(2);
                
                %Calculate distances of all particles in first frame to user
                %input points
                distances=NaN(obj.NOP,2);
                tr1=trace(find(trace(:,3)==1),:);
                for p=1:obj.NOP
                    trp=tr1(find(tr1(:,4)==p),:);
                    distances(p,1)=sqrt((trp(1,1)-x(1))^2+(trp(1,2)-y(1))^2);
                    distances(p,2)=sqrt((trp(1,1)-x(2))^2+(trp(1,2)-y(2))^2);
                end
                %Find minimal distances and corresponding particle IDs
                min1=nanmin(distances(:,1));
                min2=nanmin(distances(:,2));
                ID1=find(distances(:,1)==min1);
                ID2=find(distances(:,2)==min2);
                
                %Get trace of particle IDS
                tr_1=trace(find(trace(:,4)==ID1),:);
                tr_2=trace(find(trace(:,4)==ID2),:);
                tracenew=[tr_1(:,1:3) ones(length(tr_1(:,1))) ; tr_2(:,1:3) 2*ones(length(tr_1(:,1)))];
                posnew=[tr_1(:,1:3) ; tr_2(:,1:3)];
            else
                error('You only have one particle. That is never going to work.');
            end
            
            %Create dropletpair object as output, with the right traces and
            %the parameters from this object
            newname=[obj.foldername obj.filename '.tif'];
            params=get_parameters(obj);
            output=dropletpair2(newname,params,'n');
            output.scale=obj.scale;
            output.framerate=obj.framerate;
            output=add_pos(output,posnew);
            output=add_trace(output,tracenew);
            output.im=obj.im;
            output.NOP=2;
            
            plot_tracks(output);
        end
        function obj=contour(obj,npoints,factor)
            
            % PIV_output : n*3 matrix containing position and radius of particle
            % npoints: number of points for contour (default=20)
            % factor: percentage of radius used for contour (default=1.1)
            
            % Input: npoints: number of points to use for contour formatoin
            %       factor: multiplied by radius to get mask size (1 equals radius,
            %       larger is recommended to get a bit of overlap)
            %       PIV_output
            %           size: n*3 with n is number of frames
            %           [x,y,r] with x,y are coordinates of particle each frame (from
            %           track)
            %           r is size of particle each frame (from circlefinder)
            
            % Overrides existing contour file
            
            %l=length(PIV_output);
            %Set number of points used to build circular mask
            %npoints=20;
            %set overlap factor for mask: the percentage of radius you want the mask to
            %overlap
            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            
            if nargin==1
                npoints=20;
                factor=1.1;
            elseif nargin==2
                factor=1.1;
            end
            
            PIV_output=[obj.tr(:,1:2) ones(length(obj.tr(:,1)),1)];
            l=length(obj.tr(:,1));
            
            ymask=NaN(npoints,l);
            xmask=NaN(npoints,l);
            for i=1:l
                x=PIV_output(i,1);
                y=PIV_output(i,2);
                r=PIV_output(i,3)*factor;
                x1=NaN(1,npoints);
                y1=NaN(1,npoints);
                for j=1:npoints
                    x1(j)=x+cos(2*pi*j/npoints)*r;
                    y1(j)=y+sin(2*pi*j/npoints)*r;
                end
                ymask(:,i)=y1;
                xmask(:,i)=x1;
            end
            save('contour.mat', 'xmask', 'ymask');
            obj.contourfile=cell(2,1);
            obj.contourfile{1,1}=xmask;
            obj.contourfile{2,1}=ymask;
            
        end
        function output=pick_selection(obj,n)
            %Shows the first image of a movie and lets your pick a number
            %of droplets that are saved into a subselection based on their
            %position
            %Stores as new object of type droples.
            %WARNING: DO NOT OVERRIDE OBJECT WITH THIS FUNCTION
            
            %Check input
            if nargin<2
                message='How many particles do you want to pick: ';
                n=get_numerical_input(message);
            end            
            if isempty(obj.tr)
                obj=track_obj(obj);
            end
            if isempty(obj.im)
                error('This function requires an image. Please add initial frame first.')
            end
            
            %plot tracks
            plot_tracks(obj);
            
            %Obtain image
            A=obj.im;
            
            %show image and ask user for particle positions
            figure(100)
            colormap(gray),imagesc(A);
            disp('Please select the droplets for the selection by clicking on them');
            [x0,y0]=ginput(n);
            close 100
            
            %Obtain particle positions
            pos_1=obj.tr(find(obj.tr(:,3)==1),:);
            ID=zeros(n,1);
            count=0;
            tr_sel=[];
            %Compare particle positions to user chosen positions
            for i=1:n
                d=sqrt((pos_1(:,1)-x0(i)).^2+(pos_1(:,2)-y0(i)).^2);
                ID_temp=obj.tr(find(d==min(d)),3);
                %if not redundant, store id
                if ~any(ID==ID_temp)
                    count=count+1;
                    tr_i(1:obj.NOF,1:3)=obj.tr((ID_temp-1)*obj.NOF+1:ID_temp*obj.NOF,1:3);
                    tr_i(1:obj.NOF,4)=count;
                    tr_sel=[tr_sel ; tr_i];
                    ID(count)=ID_temp;
                end
            end
            %remove zeroes
            ID=ID(1:count);
            %remove NaN
            ID(isnan(ID))=[];
            %save selection into new object of type droplets
            fname=[obj.foldername obj.filename '.tif'];
            parameters=get_parameters(obj);
            
            output=droplets(fname,parameters,'n');
            output=change_scale(output,obj.scale);
            output=change_framerate(output,obj.framerate);
            output=add_trace(output,tr_sel);
            
        end
    end
end