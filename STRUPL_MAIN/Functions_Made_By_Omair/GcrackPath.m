function crack_path_2 = GcrackPath(Normal_stress, Number_of_Nodes, nodal_connectivity_values, sigma_t, Coordinates_Plate,Purpose)

%% Function used to calculate the nodal Stresses from a matrix with element stresses
[rows,columns] = size(Normal_stress);

Repitition_Remover = Node_Repitition_Remover(nodal_connectivity_values);

crack_path_1 = Node_Highest_Stress_Identifier(nodal_connectivity_values,Normal_stress,sigma_t,Purpose);

crack_path_2 = zeros(21,3);

%% Getting Element Numbers too

for k = 1:Number_of_Nodes
    for i = 1:rows
        for j = 1:columns
            if nodal_connectivity_values(i,j) == k
                if Repitition_Remover(i,j) == 1
                    crack_path_2(k,1) = nodal_connectivity_values(i,j);
                    crack_path_2(k,2) = crack_path_1(i,j);
                end
            end
        end
    end
end

Cracks_All = {};
node_configuration = Node_Configurator(Coordinates_Plate);

for k = 1:Number_of_Nodes %going over the nodes
    if crack_path_2(k,2) ~= 0 %chekcing if the point is overstressed
        Crack_Path_Draw = [];
        Length = 0;
        
        for l = k-3:k+3 %Looping over nearby nodes
            if l > 1 && l < Number_of_Nodes
                if crack_path_2(l,2) ~= 0 %Checking if the nearby nodes are overstressed
                    % Finding Differences and adding them up
                    Length = Length + (0.5) * sqrt((Coordinates_Plate(k,1) - Coordinates_Plate(l,1))^2 +...
                        (Coordinates_Plate(k,2) - Coordinates_Plate(l,2))^2);
                    Crack_Path_Draw(l-1,1) = l;
                    Crack_Path_Draw(l-1,2) = Coordinates_Plate(l,1);
                    Crack_Path_Draw(l-1,3) = Coordinates_Plate(l,2);



                end

            end
        end
        crack_path_2(k,3) = Length; %%assisning additions to the matrix
        [row_Cracks,column_Cracks] = size(Cracks_All);
        Cracks_All(row_Cracks + 1,1) = {Crack_Path_Draw};
    end
end


    
end