function output=correlate_functions(f,g)
%Takes to inputs functions f and g as a function of time.
%f and g need to be organized as [t q] where t is the time and q is the
%variable of interest. f and g need to have the same size. f and g need to
%have the same spacing in time.
%f and g need to be continuous in time (no gaps);
%Calculates the correlation function for a maximum time interval of 10% of
%the total measured time

%Created on 

%assert that f and g have the same size
if isempty(f) || isempty(g)
    output=NaN;
    disp('One or more of the input functions is emtpy. The correlation function does not exist');
elseif size(f)~=size(g)
    output=NaN;
    disp('The two input functions have different sizes. No correlation function can be calculated');
end

ll=length(f(:,1)); %Length of the matrix
tmax=floor(length(f(:,1))/10); %Find max time interval
tmin=1; %Min time interval is 1.
corr=NaN((tmax-tmin)+1,2); %Create empty matrix to save correlation function
corr(1,1:3)=[0 0 0]; %The correlation at t=0 is 0.

%for all time intervals
for t=tmin:tmax
    %shift f and g by the time interval
    ft=f(t+1:t:ll,2);
    f0=f(1:t:ll-t,2);
    gt=g(t+1:t:ll,2);
    g0=g(1:t:ll-t,2);
    %correlate the time-shifted functions
    corrt=(ft-f0).*(gt-g0);
    %time average
    corr(t+1,1)=t;
    corr(t+1,2)=nanmean(corrt);
    corr(t+1,3)=nanstd(corrt);
end

%store output
output=corr;

end
