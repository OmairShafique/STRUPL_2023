function Reduced_Stresses_F = Stress_Reducer(Stress_Matrix,nodes_per_element,number_of_nodes,nodal_connectivity_values,sigma_t,Coordinates_Plane)

%% Stress Reducer Taken Stresses at the element level and divides them amongst the nodes, then it takes
% stresses at the nodes form individual elements and adds them togeather
[rows,columns] = size(Stress_Matrix);

Reduced_Stresses = zeros(rows,nodes_per_element);

[rows,columns] = size(Reduced_Stresses);

for i = 1:rows
    Average = Stress_Matrix(i,1)/3;
    for j = 1:columns
        Reduced_Stresses(i,j) = Average; %forming a matrix that will provide stresses at each node 
        % from averging the stress at the element
    end
end

%% This porting utilizes Gcrack to provide stresses added at each node
Reduced_Stresses_F = Gcrack(Reduced_Stresses,number_of_nodes,nodal_connectivity_values,sigma_t,Coordinates_Plane,'Reduction');

%% This portion removes the last column so we don't have 
Reduced_Stresses_F(:,3) = [];

end