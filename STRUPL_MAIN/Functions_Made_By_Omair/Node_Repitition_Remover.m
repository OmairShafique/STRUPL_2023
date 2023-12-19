function nodal_connectivity_Repitition_Remover = Node_Repitition_Remover(nodal_connectivity_values)
%% This function removes any node repitition and provides a matrix with 1's where the nodes occur first
%and 0's where they occurs twice or 3rice

[rows,columns] = size(nodal_connectivity_values);

%% Logic Starts here
% Initiating an empty elements x nodes-per-element matrix
nodal_connectivity_Repitition_Remover = zeros(rows,columns);
index_length = rows*columns;

% Taking transpose
nodal_connectivity_Repitition_Remover = nodal_connectivity_Repitition_Remover';
nodal_connectivity_values = nodal_connectivity_values';

% Taking the first 3 constituents to be constant
nodal_connectivity_Repitition_Remover(1:3,1)  = nodal_connectivity_values(1:3,1);

% Providing a variable that will serve as a guide to the code to check if
% the value is repeated.
check_if_same = false;

% Starting from 4 and then loping all over the matrix to check if it
% repeats
for index = 4:index_length
    check_if_same = false;
    for index_inside = 1:index-1
        if check_if_same
            %if it repeated last cycle don't go ahead continue loop to next node.
            continue;
        end
        if nodal_connectivity_values(index) == nodal_connectivity_values(index_inside)
            check_if_same = true;
            % if repitition found continue and check variable
            continue;
        end
    end
    if check_if_same == false
        nodal_connectivity_Repitition_Remover(index) = nodal_connectivity_values(index);
        % else assign the node value at that location.
    end
    
end

nodal_connectivity_Repitition_Remover = nodal_connectivity_Repitition_Remover';

for i = 1:rows 
    for j = 1:columns
        if nodal_connectivity_Repitition_Remover(i,j) ~= 0
            nodal_connectivity_Repitition_Remover(i,j) = 1;
            %change node number to 1's and 0's
        end
    end
end

end

