classdef chains < particles
    %Chains class is a sublcass of particles unique for having one
    %dimensional aggregates of particles in the images
    %It has no extra properties with respect to particles but overrides
    %some methods and has some new methods
    %
    %Make sure this code is in the same folder as particles.m and all
    %helper codes
    %
    %Extra methods
    % - obj=chains(fname,parameters,identity)       Inintializes object of type chains 
    % - output=get_parameters(obj)                  overrides particles function 
    % - obj=read_parameters(obj,parameters)         overrides particles function 
    % - split_to_single_chain(obj)                  finds individual 1D aggregates and creates single_chain objects for each 
    %
    %Created on 19-07-17 by Pepijn Moerman
    %Last modified: 26-02-18
    
    properties
        %no new properties in this object
    end
    methods
        function obj=chains(fname,parameters,identity)
            %Chains needs to be an object of type particles
            %Identity indicates wether or not analysis should be started
            %'y' is start, anything else is not.
            if nargin==1
                parameters=[];
            
            elseif nargin<1
                parameters=[];
                fname=[];
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
                answer2='n';
                while answer=='y'
                    if answer1=='n'
                        obj=find_pos(obj);
                        message='Are you happy with the particle finding? Yes (y) or no (n): ';
                        answer1=get_textual_input(message,acceptables);
                    end
                    if answer1=='y' && answer2~='y'
                        proceed=1;
                        try obj=track_obj(obj);
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
                                answer='y';
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
            parameters{5,1}='cutoffvalue';
            parameters{5,2}=obj.cutoff;
        end
        function obj=read_parameters(obj,parameters)
            %parameters needs to be a 4x2 cell created using the
            %get_parameters function of this object
            obj.framerate=parameters{1,2};
            obj.scale=parameters{2,2};
            obj.posparam=parameters{3,2};
            obj.trparam=parameters{4,2};
            obj.cutoff=parameters{5,2};
        end
        function obj=find_neighbours(obj,cutoffvalue)
            %creates matrix num_neighs(indicating how many neighbours each particle has)
            %creates matrix neighbours (indicating which particles are the
            %neighbours)
            %neighbours are defined using a cutoff value
            %Assumes that the number of neighbours is a fixed value for one
            %object and does not change from frame to frame
            
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
                neighbourstest=find(av_dist(i,:)<cutoffvalue);
                neighbours_temp=neighbourstest(find(neighbourstest~=i));
                obj.neighbours{i}=neighbours_temp;
                num_neighs_temp(i,1)=length(neighbourstest)-1;
            end
            obj.num_neighs=num_neighs_temp;
            obj.cutoff=cutoffvalue;
        end
        function singles=split_to_single_chain(obj)
            %uses neighours to split a set of multiple chains into single
            %chains.
            %outputs a cell with length 'number of chains' and in each cell
            %element an object of type 'single_chain'
            %ignores all single particles, doublets and branched objects;
            
            %Checks if neighbours is empty, otherwise fills neighbours
            if isempty(obj.neighbours)
                obj=find_neighbours(obj);
            end
            
            %creates a temporary copy of neighbours
            neighbours_temp=obj.neighbours;
            %creates a temporary number of neighbours
            nn_temp=obj.num_neighs;
            %create an empty cell to save individual chains
            chains=cell(obj.NOP,1);
            %count keeps track of number of chains
            count=0;
            %get parameters and folder for saving into new object later
            fname=[obj.foldername obj.filename '.tif'];
            parameters=get_parameters(obj);
            
            %loop through all particles
            for i=1:obj.NOP
                nogo=1;
                %check if particle was previously incorporatded into a chain
                for j=1:i-1
                    if any(chains{j}==i)
                        nogo=0;
                    end
                end
                %if particle is not previously incorporated
                if nogo==1
                    %chains temp is 1D matrix containing ID's of all particles in that chain
                    chains_temp=findchain_class(neighbours_temp,nn_temp,i);
                    %if chain is longer than 2, keep, otherwise discard
                    if length(chains_temp)>1 %! If you are interested in dimers or monomers also, change this line to if ch_te not equal NaN
                        %increase number of chains found
                        count=count+1;
                        %save chains_temp into a cell
                        chains{count}=chains_temp;
                    end
                end
            end
            
            %create a cell of single_chains. This is the output variable
            singles=cell(count,1);
            %for all found chains
            for i=1:count
                %create a single chain object with name and parameters
                output=single_chain(fname,parameters,'n');
                %save the length of this chain
                ll=length(chains{i});
                %create a temporary list, chains temp again
                chains_temp=chains{i};
                %create an empty matrix for the trace of all the particles in this chain
                tr_temp=[];
                %create an empty cell of neighbours with length of chain
                neighbours_j=cell(1,ll);
                %create an empty cell of the number of neighbours with length of chain.
                numneighs_j=NaN(ll,1);
                %for all particles in chain
                for j=1:ll
                    %convert old ID into new ID with increasing value for particles in one chain.
                    ID_old=chains_temp(j);
                    ID_new=j;
                    tr_j=[obj.tr(find(obj.tr(:,4)==ID_old),1:3) ones(obj.NOF,1)*ID_new];
                    tr_temp=[tr_temp ; tr_j];
                    nn_j=neighbours_temp{ID_old};
                    %Relabel neighbours_j
                    if j==1 || j==ll
                        numneighs_j(j,1)=1;
                        nn_j(1)=find(chains_temp==nn_j(1));
                    else
                        numneighs_j(j,1)=2;
                        nn_j(1)=find(chains_temp==nn_j(1));
                        nn_j(2)=find(chains_temp==nn_j(2));
                    end
                    neighbours_j{j}=nn_j;
                end
                %Run the single_chain machinery over it
                output=add_trace(output,tr_temp);
                output.NOP=ll;
                output.neighbours=neighbours_j;
                output.num_neighs=numneighs_j;
                output.chainlength=ll;
                output=make_dist_matrix(output);
                output.cutoff=obj.cutoff;
                %output=find_angles(output);
                output=find_center_of_mass(output);
                %output=calc_diff_tensor(output,1);
                %output=find_endvector(output);
                %output=find_angle_MSD(output);
                singles{i}=output;
            end %end count
            
        end
        function singles=split_to_dimer(obj)
            %uses neighours to split a set of multiple chains into dimers
            %outputs a cell with length 'number of chains' and in each cell
            %element an object of type 'dimer'
            %ignores all clusters of other sizes
            
            %Checks if neighbours is empty, otherwise fills neighbours
            if isempty(obj.neighbours)
                obj=find_neighbours(obj);
            end
            
            %creates a temporary copy of neighbours
            neighbours_temp=obj.neighbours;
            %creates a temporary number of neighbours
            nn_temp=obj.num_neighs;
            %create an empty cell to save individual chains
            chains=cell(obj.NOP,1);
            %count keeps track of number of chains
            count=0;
            %get parameters and folder for saving into new object later
            fname=[obj.foldername obj.filename '.tif'];
            parameters=get_parameters(obj);
            
            %loop through all particles
            for i=1:obj.NOP
                nogo=1;
                %check if particle was previously incorporatded into a chain
                for j=1:i-1
                    if any(chains{j}==i)
                        nogo=0;
                    end
                end
                %if particle is not previously incorporated
                if nogo==1
                    %chains temp is 1D matrix containing ID's of all particles in that chain
                    chains_temp=findchain_class(neighbours_temp,nn_temp,i);
                    %if chain is longer than 2, keep, otherwise discard
                    if ~isnan(chains_temp) %! If you are interested in dimers or monomers also, change this line to if ch_te not equal NaN
                        if length(chains_temp)==2
                            %increase number of chains found
                            count=count+1;
                            %save chains_temp into a cell
                            chains{count}=chains_temp;
                        end
                    end
                end
            end
            
            %create a cell of single_chains. This is the output variable
            singles=cell(count,1);
            %for all found chains
            for i=1:count
                %create a single chain object with name and parameters
                output=dimer(fname,parameters,'n');
                %save the length of this chain
                ll=length(chains{i});
                %create a temporary list, chains temp again
                chains_temp=chains{i};
                %create an empty matrix for the trace of all the particles in this chain
                tr_temp=[];
                %create an empty cell of neighbours with length of chain
                neighbours_j=cell(1,ll);
                %create an empty cell of the number of neighbours with length of chain.
                numneighs_j=NaN(ll,1);
                %for all particles in chain
                for j=1:ll
                    %convert old ID into new ID with increasing value for particles in one chain.
                    ID_old=chains_temp(j);
                    ID_new=j;
                    tr_j=[obj.tr(find(obj.tr(:,4)==ID_old),1:3) ones(obj.NOF,1)*ID_new];
                    tr_temp=[tr_temp ; tr_j];
                    nn_j=neighbours_temp{ID_old};
                    %Relabel neighbours_j
                    if j==1 || j==ll
                        numneighs_j(j,1)=1;
                        nn_j(1)=find(chains_temp==nn_j(1));
                    else
                        numneighs_j(j,1)=2;
                        nn_j(1)=find(chains_temp==nn_j(1));
                        nn_j(2)=find(chains_temp==nn_j(2));
                    end
                    neighbours_j{j}=nn_j;
                end
                %Run the dimer machinery over it
                output=add_trace(output,tr_temp);
                output.NOP=ll;
                output.neighbours=neighbours_j;
                output.num_neighs=numneighs_j;
                output.chainlength=ll;
                output=make_dist_matrix(output);
                output.cutoff=obj.cutoff;
                output=find_center_of_mass(output);
                %output=calc_diff_tensor(output,1);
                %output=find_endvector(output);
                singles{i}=output;
            end %end count
            
        end
    end
end