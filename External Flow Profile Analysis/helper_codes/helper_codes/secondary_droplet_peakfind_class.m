function output=secondary_droplet_peakfind_class(A,fname,k,find_parameters,dist_filt_check)

% Unpack information in find_parameters
diameter = find_parameters{1};
lnoise = find_parameters{2};
r_lower = find_parameters{3};
r_upper = find_parameters{4};
sensitivity = find_parameters{5};
invert_check=find_parameters{6};
pk_old=find_parameters{7};

% Process data with values from initial image
if invert_check=='y'
    A=255-A;
end
Ai=bpass(A,lnoise,diameter);
[centers, radii] = imfindcircles(Ai,[round(r_lower),round(r_upper)],'Sensitivity', sensitivity);

% Show the image with the newly found peaks
if isempty(centers)
    pkl=0;
else
    pkl=length(centers(:,1));
end
figure
hold on
colormap('gray'),imagesc(Ai);
if ~isempty(centers)
    viscircles(centers(:,1:2),radii(:,1));
    pk=centers(:,1:2);
else
    pk=NaN;
end
hold off

    

% Make sure we find the right particles
message='These are the particles found in the second frame using the parameters you optimized for the first frame. Did it find the particles correctly? (y/n): ';
acceptables=['y' 'n'];
try_again_check = get_textual_input(message,acceptables);

if try_again_check=='n'
    disp('Sadly the parameters you chose initially do not seem to work for the next frame. This does not bode well for the rest of the data.')
    disp('To improve the situation we are now going to find the parameters with which we can find the particles in this frame')
    % If image was inverted for inital optimization, invert back so we can
    % start second optimization with a clean slate
    A = imread(fname, k);
    new_parameters=initial_droplet_peakfind_class(A);
    pk=new_parameters{7};
elseif try_again_check=='y'
    disp('Hurray, the parameters you chose worked for two consecutive frames boding well for the rest of your data.')
end

if dist_filt_check==1
% Calculate the shortest distances from frame 1 to frame 2 and use those
% distances to obtain an estimate of the maximum allowed distance

if length(pk(:,1))==length(pk_old(:,1))
    dist=NaN(length(pk(:,1)),length(pk_old(:,1)));
    min_dist=NaN(length(pk(:,1)),1);
    for i=1:length(pk(:,1))
        for j=1:length(pk_old(:,1))
            dist(i,j)=sqrt((pk(i,1)-pk_old(j,1))^2+(pk(i,2)-pk_old(j,2))^2);
        end
        min_dist(i,1)=min(dist(i,:));
    end
else
    error('There are more particles found in your second frame than your first frame. The code cannot deal with that yet')
end

distmin=3*max(min_dist(:,1));
disp(['The maximum displacement of your particles was ' num2str(distmin) '.'])
disp('This value is saved and used to throw away particles found more than 3 times this distance from a previous position');
else
    distmin=5;
end

%Package the new parameter sets
new_parameters{8}=distmin;
find_parameters{8}=distmin;
sum_parameters{1}=find_parameters;
if try_again_check=='n'
    sum_parameters{2}=new_parameters;
    disp('ID2')
    new_parameters
end


output={distmin pk sum_parameters};