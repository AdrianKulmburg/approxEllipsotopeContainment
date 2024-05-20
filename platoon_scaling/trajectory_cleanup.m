function y = trajectory_cleanup(y, tol)
% When plotting the trajectories and transforming them into .tex files,
% there may be too many points, which overloads latex. Therefore, we need
% to selectively remove some points that are too close to one another, as
% in this case they are pretty redundant.

current_index = 2;
while current_index <= size(y,2)
    y1_previous = y(1,current_index-1);
    y2_previous = y(2,current_index-1);
    
    y1 = y(1,current_index);
    y2 = y(2,current_index);
    
    if (abs(y1-y1_previous) > tol) || (abs(y2-y2_previous) > tol)
        current_index = current_index + 1;
    else
        y(:,current_index) = [];
    end


end