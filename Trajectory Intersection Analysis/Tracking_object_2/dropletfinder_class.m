function output=dropletfinder_class(fname,parameters)%fname contains name of file with tifstack
%Finds droplet centers in a tiffstack using hough transform
%

% Modification History:
% Based on particle_finder_class
% and track_and size
%   Written by Pepijn 11-2016
%   Modified by Pepijn 01-09-2017


%% Initialize

% Make sure the number of input variables is correct
if nargin<2
    initialize=1;
else
    initialize=0;
    sum_parameters=parameters{1};
    answers=parameters{2};
    imageshow=answers{1};
    dist_filt_check=answers{2};
    num_filt_check=answers{3};
end

% Set warning for circlefind radius off.
MSGID1='images:imfindcircles:warnForSmallRadius';
warning('off', MSGID1)
% Set warning non integer index off
MSGID2='MATLAB:colon:nonIntegerIndex';
warning('off', MSGID2)

%assert good value for num_images
info = imfinfo(fname);
num_images = numel(info);

%Define matrix
pos=[];

%% Loop through stack and find droplets

count=0;
disp(['Finding particles in ' num2str(num_images) ' frames.'])
%for all images in stack
for k = 1:num_images
    %read image
    
    A = imread(fname, k);
    %print image number
    if ~mod(k/100,1)
        disp(['Working on frame: ' num2str(k) '.'])
    end
    
    %Assert image is 8-bit
    if length(size(A))>2
        disp('Image is not of type 8-bit. Converting to 8bit first.')
        %Convert to grayscale
        A = rgb2gray(A);
    end
    
    %Stretch contrast to go from 1-255;
    A = uint8((double(A)-double(min(A(:))))/(double(max(A(:)))-double(min(A(:))))*255);
       
    %For the first image of the first movie
    if count==0 && initialize==1
        count=count+1;
        disp('Hello user. We are going to do some parameter optimization before starting particle tracking')
        find_parameters=initial_droplet_peakfind_class(A);
        pk=find_parameters{7};
        
        
        % Ask for filters
        message= 'Do you want to apply a distance_filter? This will ensure you only find particles close to particles in the previous frame. Yes (y) or no (n): ';
        acceptables=['y' 'n'];
        dist_filt_check_text=get_textual_input(message,acceptables);
        if dist_filt_check_text=='y'
            dist_filt_check=1;
        else
            dist_filt_check=0;
        end

        message= 'Is the number of particles you wish to find the same in each frame of each movie?. Yes (y) or no (n): ';
        acceptables=['y' 'n'];
        num_filt_check_text=get_textual_input(message,acceptables);
        if num_filt_check_text=='y'
            num_filt_check=1;
        else
            num_filt_check=0;
        end

        %For the second image of the first movie: check and get distance_limit
        
    elseif count==1 && initialize==1
        count=count+1;
        aa=secondary_droplet_peakfind_class(A,fname,k,find_parameters,dist_filt_check);
        dist_limit=3*aa{1};
        pk=aa{2};
        sum_parameters=aa{3};
        
        %For all the other images
    else
        
        output=dropletoptimizer_class(A,sum_parameters,k,pos,dist_filt_check,num_filt_check);
        pk=output{1};
        sum_parameters=output{2};
        if  count==2
            count=count+1;
            message= 'Do you want to show images with found particles? If so, every how many frames? Answer 0 if no frames: ';
            imageshow=get_numerical_input(message);
        end
        %Plot image with found particles if imageshow=1
        if imageshow~=0
            if ~mod(k/imageshow,1) && ~isempty(pk)
                find_parameters=sum_parameters{1};
                r_lower=find_parameters{3};
                pkl=length(pk(:,1));
                radiplot=zeros(pkl,1)+r_lower;
                figure(k)
                hold on
                colormap('gray'),imagesc(A);
                viscircles(pk(:,1:2),radiplot);
                hold off
            end
        end
        
        
    end
    
    %Polish the matrix containing particle coordinates
%     if isempty(pk)
%         pk(1,1:2)=0;
%     else
%         pk=pk(:,1:2);
%     end
    try pk=pk(:,1:2);
    catch
        pk(1,1:2)=0
    end

    pk(:,3)=k;
    pos=[pos ; pk];
    
    
end

% Remove particles found at the origin (these can be an artefact of tracking)
pos(find(pos(:,1)==0),:)=[];

answers={imageshow,dist_filt_check,num_filt_check};
parameters={sum_parameters answers};

output={pos parameters};
end
