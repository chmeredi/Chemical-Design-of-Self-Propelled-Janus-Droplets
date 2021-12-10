function output = findchain_class(neighbours,num_neighs,i)

%Function Name: findchain
%Goal: finds the length of a 1D aggregate of particles attached to a
%   intial particle
%Expected input:
%   neighbours: a 1xn cell with n is the number of particles. Each cell
%   element is a matrix containing which particles are the neighbours
%   num_neighs is a 1xn matrix with each element containing the number
%   of neighbours
%   i: initial particle index: integer
%Output: a 1xm matrix containing the ID's of all particles in the
%chain, where m is the length of the chain
%        Output is '0' if aggregate is not 1D or has branches
%Created by Pepijn Moerman on 07-07-16
%Last modified on: 15-01-18 by Pepijn Moerman

%create an empty list of ID's that we will fill with the particles in the chain.
output=[];

%if particle i is the start of a chain
if num_neighs(i)==1
    output=[output i];
    %find the neighbour of particle i;
    i=neighbours{i};
    %For as long as the chain continues
    while 1==1
        %if the particle has only one neighbour, you've found the end
        if num_neighs(i)==1
            output=[output i];
            % break while loop
            break
        %if there is more than 2 neighbours, this is not a linear cluster, abort
        elseif num_neighs(i)~=2
            %create empty output
            output=[];
            %and break out of loop
            break
        %otherwise there is a chain continuation, then redo loop
        else
            %Add the ID to the list of IDs of this chain
            output=[output i];
            %find the neighbours of particle i;
            neighs_temp=neighbours{i};
            %remove neighbours equal to previous particles
            neighs_temp(find(neighs_temp==output(end-1)))=[];
            %Make the new particle i the neighbour
            i=neighs_temp;
        end
    end %end while loop
end %end check initial particle has 1 neighbour

end %end function




