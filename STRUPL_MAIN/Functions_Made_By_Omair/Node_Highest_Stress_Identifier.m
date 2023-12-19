function node_highest_stress = Node_Highest_Stress_Identifier(nodal_connectivity_values,Normal_stress,sigma_t,Purpose)

Repitition_Remover = Node_Repitition_Remover(nodal_connectivity_values);
Repitition_Remover = Repitition_Remover';

[rows,columns] = size(nodal_connectivity_values);

node_highest_stress = zeros(rows,columns);
index_length = rows*columns;


node_highest_stress = node_highest_stress';
nodal_connectivity_values = nodal_connectivity_values';
Normal_stress = Normal_stress';




for index = 1:index_length
    max_value = 0;
    % check_if_same = false;
    for index_inside = 1:index_length
        % if check_if_same
        %     continue;
        % end 
            if nodal_connectivity_values(index) == nodal_connectivity_values(index_inside)
                if abs(Normal_stress(index_inside)) > max_value
                    switch Purpose
                        case 'Reduction'
                            max_value = max_value + Normal_stress(index_inside);
                        case 'Crack_path'
                            max_value = Normal_stress(index_inside);
                    end
                    
                end
            end
    end
    if Repitition_Remover(index) == 1
        if abs(max_value) > abs(sigma_t)
            node_highest_stress(index) = max_value;
        end
    end
end

    
    % if check_if_same == false

node_highest_stress = node_highest_stress';

% for i = 1:rows 
%     for j = 1:columns
%         if node_highest_stress(i,j) ~= 0
%             node_highest_stress(i,j) = 1;
%         end
%     end
% end

end

