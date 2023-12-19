function fg_gravity_Final_end = fg_matrix_calculator(number_of_elements,nodal_connectivity_values,number_of_nodes,nf)

%% Finding a matrix that will provide the number of times a node is repeated and work with fg_matrix to 
% find the gravity loading matrix in tune witht the Global_Force_Vector so
% these two can be added.

[rows,columns] = size(nodal_connectivity_values);
index_length = rows*columns;



% node_repitition = zeros(number_of_nodes,1);
% 
% for i = 1:number_of_nodes
%     times_repeated = 0;
%     for index = 1:index_length
%         if i == nodal_connectivity_values(index)
%             times_repeated = times_repeated + 1;
%         end
%     end
%     node_repitition(i) = times_repeated;
% end
% 
% node_repitition_modified = nf .* node_repitition;
% 
% [rows,columns] = size(node_repitition_modified);
% index_length = rows*columns;
% 
% node_repitition_final = zeros(index_length,1);
% 
% index = 1;
% for i = 1:rows
%     for j = 1:columns
%         node_repitition_final(index) = node_repitition_modified(i,j);
%         index = index + 1;
%     end
% end


gamma = evalin('base','weigth_density'); 
thickness_of_plate = evalin('base','thickness_of_plate');
number_of_dof_per_node = evalin('base','nodof');
Element_Type = evalin('base','Element_Type');
number_of_nodes_per_element = evalin('base','number_of_nodes_per_element');
nodal_coordinate_values = evalin('base','nodal_coordinate_values');
Number_of_Nodes = length(nodal_coordinate_values); % Infered Number of Nodes in the System
fg_gravity = zeros(number_of_elements * number_of_nodes_per_element * number_of_dof_per_node,1);

if Element_Type == 3
    number_of_nodes_per_element = 3;
else
    number_of_nodes_per_element = length(nodal_connectivity_values(1,:)); %Number of Nodes per element
end

Gravity_Load = [0 ,-gamma]';

current_row = 1;
for iel = 1:number_of_elements

        [bee,fun_3,g,A,d_3] = elem_T3(iel);

        fg_gravity(current_row: current_row + number_of_dof_per_node*number_of_nodes_per_element - 1) = ...
            (fun_3 * Gravity_Load * d_3) * thickness_of_plate * (-1/3);

        current_row = current_row + number_of_dof_per_node * number_of_nodes_per_element;
end

[rows,columns] = size(fg_gravity);

fg_gravity_x2 = zeros(rows/2,3);

index = 1;

for i = 1:rows/2
    for j = 2:3
        fg_gravity_x2(i,j) = fg_gravity(index);
        index = index + 1;
    end
end

[rows,columns] = size(nodal_connectivity_values);
index = 1;

for i = 1:rows
    for j = 1:columns
        fg_gravity_x2(index,1) = nodal_connectivity_values(i,j);
        index = index + 1;
    end
end

fg_gravity_Final = zeros(Number_of_Nodes,3);
[rows,columns] = size(fg_gravity_x2);
index = 1;

for i = 1:Number_of_Nodes
    for j = 1:rows
        if fg_gravity_x2(j,1) == i
            fg_gravity_Final(i,1) = i;
            fg_gravity_Final(i,2:3) = fg_gravity_Final(i,2:3) + fg_gravity_x2(j,2:3);
        end
    end
end

[rows,columns] = size(fg_gravity_Final);
index = 1;
fg_gravity_Final_end = zeros(rows*(columns-1),1);

for i = 1:rows
    for j = 2:columns
        fg_gravity_Final_end(index) = fg_gravity_Final(i,j);
        index = index + 1;
    end
end


%Adding ____ load to fg matrix


end