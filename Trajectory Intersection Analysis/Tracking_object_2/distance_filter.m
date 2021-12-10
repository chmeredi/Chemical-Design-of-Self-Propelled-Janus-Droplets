function output=distance_filter(pos_new,pos_ref,dist_limit)

% Takes 2 lists of positions of format nx2 [xi yi] and a value distance limit
% Removes positions from pos_new that are not within a distance dist_limit
% from any position is pos_ref


new_positions=[];
count=0;
for i=1:length(pos_new(:,1))
    x_new=pos_new(i,1);
    y_new=pos_new(i,2);
    distance=NaN(length(pos_ref(:,1)),1);
    for j=1:length(pos_ref(:,1))
        x_ref=pos_ref(j,1);
        y_ref=pos_ref(j,2);
        distance(j,1)=sqrt((x_ref-x_new)^2+(y_ref-y_new)^2);
    end
    mindist=min(distance(:,1));
    if mindist<dist_limit
        count=count+1;
        new_positions(count,1:2)=pos_new(i,1:2);
    end    
end

output=new_positions;