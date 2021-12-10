function output=initial_particle_peakfind_class(A)

% Ask user for input to track the particles in the intial frame correctly:
% Takes input: images A
% Produces output: values for tracking as a list
% Output = [diameter lnoise filtervalue featuresize threshold, pk]

% Uses helper codes: - get_numerical_input.m
%                    - get_textual_input.m
%                    - pkfind.m (by Crocker and Grier)

% Modification History:
%   - written by Pepijn Moerman on 10-03-2017


% Get user input for diameter
message= 'What is approximately the diameter of your particle? Answer with a numerical value: ' ;
diameter= get_numerical_input(message);

% Possibly invert image
message= 'Your image looks like this. Note that this tracking software finds light spots in a dark background. Do you want to invert the contrast? ' ;
acceptables=['y' 'n'];
figure
hold on
colormap('gray'),imagesc(A);
hold off
invert_check = get_textual_input(message,acceptables);
if invert_check=='y'
    A=255-A;
end

accept_total='n'; 
while accept_total~='y'
% Get user input for bpass noise level and apply bpass
disp('Default noise level for bandpass is 1. This is the image after running the bandpass.')
accept='n';
lnoise=1;
while accept~='y' && accept_total~='g'
    Ai = bpass(A,lnoise,diameter);
    figure
    hold on
    colormap('gray'),imagesc(Ai);
    hold off
    message='Are you happy with this processing? If so: press "y". If not: press "d" to change diameter or "b" to change bandpassfilter: ';
    acceptables=['y' 'b' 'd'];
    accept=get_textual_input(message,acceptables);
    if accept=='d'
        message= 'What new value for diameter do you want to use? ' ;
        diameter= get_numerical_input(message);
    elseif accept=='b'
        message= 'What new value for bandpassfilter do you want to use? ' ;
        lnoise= get_numerical_input(message);
    end
end

% Take processed image and apply gaussian filter
disp('Next we use a gaussian filter to further impove the image. Again you need to give some input')
accept='n';
filtervalue=diameter;
while accept~='y'
    disp('The default value for the gaussian filter is the particle diameter')
    Aj=imgaussfilt(Ai,filtervalue);
    figure
    hold on
    colormap('gray'),imagesc(Aj);
    hold off
    message='Are you happy with this processing? If so: press "y". If not: press "n": ';
    acceptables=['y' 'n'];
    accept=get_textual_input(message,acceptables);
    if accept=='n'
        message= 'What new value for the guassian filter do you want to use? ' ;
        filtervalue= get_numerical_input(message);
    end
end

disp('We are done with particle preprocessing now. It is time to look if we can find the particles.')
disp('Default values for peakfind are threshold=1 and the feature size=diameter')

% Take processed image and apply peakfind
threshold=1;
feature_size=diameter;
accept='n';
while accept ~= 'y'
    pk = pkfnd(Aj,threshold,feature_size);
    if ~isempty(pk)
        pk = cntrd(Ai,pk,2*feature_size+1);
    end

    if isempty(pk)
        pkl=0;
    else
        pkl=length(pk(:,1));
    end
    radiplot=zeros(pkl,1)+feature_size/2;
    figure
    hold on
    colormap('gray'),imagesc(Aj);
    if ~isempty(pk)
    viscircles(pk(:,1:2),radiplot);
    end
    hold off
    message='Did we find all your particles and nothing else? If so: press "y". If not: press "f" to change feature_size, press "t" to change threshold: ';
    acceptables=['y' 'f' 't'];
    accept=get_textual_input(message,acceptables);
    if accept=='f'
        message= 'What new value for feature size do you want to use? ' ;
        feature_size= get_numerical_input(message);
    elseif accept=='t'
        message= 'What new value for threshold do you want to use? ' ;
        threshold= get_numerical_input(message);
    end

end

message='Are you done with initialization? If so answer "y". If you want to go back to the bandpass, answer "b". If you want to go back to the gaussian filter, answer "g". ';
acceptables=['y' 'b' 'g'];
accept_total=get_textual_input(message,acceptables);
end

%Display parameters used
disp('You are done with tracking initialization. The following values are saved for the tracing prodecure: ')
disp(['diameter = ' num2str(diameter) '.'])
disp(['band pass noise level = ' num2str(lnoise) '.'])
disp(['gaussian filter size = ' num2str(filtervalue) '.'])
disp(['peakfind feature size = ' num2str(feature_size) '.'])
disp(['peakfind threshold = ' num2str(threshold) '.'])

close all;

output = {diameter lnoise filtervalue feature_size threshold,invert_check, pk};