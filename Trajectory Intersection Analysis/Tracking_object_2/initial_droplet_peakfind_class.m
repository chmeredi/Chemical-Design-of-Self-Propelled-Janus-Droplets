function output=initial_droplet_peakfind_class(A)

% Ask user for input to track the particles in the intial frame correctly:
% Takes input: images A
% Produces output: values for tracking as a list
% Output = [diameter lnoise filtervalue featuresize threshold, pk]

% Uses helper codes: - get_numerical_input.m (Pepijn)
%                    - get_textual_input.m (Pepijn)
%                    - imfindcircle.m (Matlab)

% Modification History:
%   - written by Pepijn Moerman on 10-03-2017
%   - modified by Pepijn Moerman on 04-09-2017


% Get user input for diameter
message= 'What is approximately the diameter of your particle? Answer with a numerical value: ' ;
diameter_actual= get_numerical_input(message);
diameter=2;

% Possibly invert image
message= 'Your image looks like this. Note that this tracking software finds light circles in a dark background. Do you want to invert the contrast? ' ;
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

disp('We are done with particle preprocessing now. It is time to look if we can find the particles.')
disp('Default values for imfindcircle are sensitivity=0.8 and the radius_lower=0.4*diameter, radius_upper=0.6*diameter')

% Take processed image and apply peakfind
sensitivity=0.9;
r_lower=0.4*diameter_actual;
r_upper=0.6*diameter_actual;
if round(r_lower<1)
    r_lower=1;
end
if round(r_upper) < 1 || round(r_upper)<round(r_lower)
    r_upper=r_lower+1;
end

accept='n';
centers=[];
radii=[];
while accept ~= 'y'
    [centers,radii] = imfindcircles(Ai,[round(r_lower),round(r_upper)],'Sensitivity',sensitivity);
    
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
    end
    hold off
    message='Did we find all your particles and nothing else? If so: press "y". If not: press "s" to change sensitivity, press "l" to change lower bound or press "u" to change upper bound: ';
    acceptables=['y' 's' 'l' 'u'];
    accept=get_textual_input(message,acceptables);
    if accept=='s'
        message= 'What new value for sensitivity do you want to use? High is sensitive. Value between 0 and 1: ' ;
        sensitivity= get_numerical_input(message);
        if sensitivity>1
            sensitivity=1;
        elseif sensitivity<0
            sensitivity=0;
        end
    elseif accept=='l'
        message= 'What new value for the lower bound of radii you want to search for? ' ;
        r_lower= get_numerical_input(message);
        if round(r_upper)<round(r_lower)
            r_lower=r_upper-1;
            if r_lower<1
                r_lower=1;
                r_upper=r_lower+1;
            end
            disp(['Your value for r_lower was lower than r_upper. We used ' num2str(round(r_lower)) ' instead.'])
        end
    elseif accept=='u'
        message= 'What new value for the upper bound of radii you want to search for? ' ;
        r_upper= get_numerical_input(message);
        if round(r_upper)<round(r_lower)
            r_upper=r_lower+1;
            disp(['Your value for r_upper was lower than r_lower. We used ' num2str(round(r_upper)) ' instead.'])
        end
    end

end

pk=centers(:,1:2);

message='Are you done with initialization? If so answer "y". If you want to go back to the bandpass, answer "b". ';
acceptables=['y' 'b'];
accept_total=get_textual_input(message,acceptables);
end

%Display parameters used
disp('You are done with tracking initialization. The following values are saved for the tracing prodecure: ')
disp(['diameter = ' num2str(diameter) '.'])
disp(['band pass noise level = ' num2str(lnoise) '.'])
disp(['lower radius bound = ' num2str(r_lower) '.'])
disp(['upper radius bound = ' num2str(r_upper) '.'])
disp(['sensitivity = ' num2str(sensitivity) '.'])

close all;

output = {diameter lnoise r_lower r_upper sensitivity invert_check pk};