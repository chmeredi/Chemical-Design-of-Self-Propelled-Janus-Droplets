function output=peakfindoptimizer_class(A,sum_parameters,k,pos,dist_filt_check,num_filt_check)

n=1;
check=0;
while check~=1
%Get the nth parameter set out of all parametersets
find_parameters=sum_parameters{n};

% Unpack information in find_parameter
if isempty(find_parameters) && n<length(sum_parameters)
    n=n+1;

elseif ~isempty(find_parameters)
    diameter = find_parameters{1};
    lnoise=find_parameters{2};
    filtervalue = 5*find_parameters{3};
    filtervalue_try=find_parameters{3};
    peakfindvalue = 5*find_parameters{4};
    peakfindvalue_try=find_parameters{4};
    threshold_try=find_parameters{5};
    threshold_max = 10 * find_parameters{5};
    invert_check=find_parameters{6};
    pk1=find_parameters{7};
    dist_limit=find_parameters{8};
    
    if ~isempty(pk1)
        n_desired=length(pk1(:,1));
    else
        n_desired=0;
        disp('No biggie, but you found no particles in the first frame. That is highly unlikely. Just FYI.')
    end

    

% Process data with values from initial image
if invert_check=='y'
    A=255-A;
end
Ai=bpass(A,lnoise,diameter);
Ai=imgaussfilt(Ai,filtervalue_try);
pk = pkfnd(Ai,threshold_try,peakfindvalue_try);
if ~isempty(pk)
    pk = cntrd(Ai,pk,2*peakfindvalue_try+1);
end

% If there is a distance constraint in the tracking
if dist_filt_check==1 && k>1
    % Apply a distance filter
    check=0;
    while check~=1
        pos_prev=pos(find(pos(:,3)==k-1),1:2);
        if ~isempty(pos_prev)
            pk_ref=pos_prev;
        else
            pk_ref=pk1;
        end
        
        pk_new=distance_filter(pk,pk_ref,dist_limit);
        
        %Make sure dist limit is not so stringent that all particles are lost
        if ~isempty(pk_new)
            check=1;
            pk=pk_new;
        elseif dist_limit>50
            error('No particles were found in this frame at all');
        else
            dist_limit=dist_limit+1;
            
        end
    end
end

% If the number of particles remains the same throughout the movie
check=1;
if num_filt_check==1
    if ~isempty(pk)
        nparticles=length(pk(:,1));
    else
        nparticles=0;
    end
    % If the number of particles is not zero and larger than n_desired
    if nparticles>n_desired && k>1
        % Keep only the particles closest to the ones in the previous frame
        pos_prev=pos(find(pos(:,3)==k-1),1:2);
        if ~isempty(pos_prev)
            pk_ref=pos_prev;
        else
            pk_ref=pk1;
        end
        pk_new=NaN(length(pk_ref),2);
        distance=NaN(length(pk(:,1)),length(pk_ref(:,1)));
        for i=1:length(pk_ref(:,1))
            distance(:,i)=sqrt((pk(:,1)-pk_ref(i,1)).^2+(pk(:,2)-pk_ref(i,2)).^2);
            [V,I]=min(distance(:,i));
            pk_new(i,1:2)=pk(I,1:2);
        end
        pk=pk_new;
    elseif nparticles<n_desired && n<length(sum_parameters)
        check=0;
        n=n+1;
        % In principle there is room for a loop through threshold to ensure
        % the right number of particles n_desired is found if the contrast
        % changes
    end
else
    if isempty(pk)&& n<length(sum_parameters)
        check=0;
    end
end

end

end%End while loop

output={pk sum_parameters};
end