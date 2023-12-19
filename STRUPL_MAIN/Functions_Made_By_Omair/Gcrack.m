function crack_path_2 = Gcrack(Normal_stress, Number_of_Nodes, nodal_connectivity_values, sigma_t, Coordinates_Plate,Purpose)

%% Function used to calculate the nodal Stresses from a matrix with element stresses
[rows,columns] = size(Normal_stress);

Repitition_Remover = Node_Repitition_Remover(nodal_connectivity_values);
% %Change to checking all nodes repittions to see if any one
% %exceeds sigma_t
%
%
%
% for i = 1:rows
%     for j = 1:columns
%         if Repitition_Remover(i,j) == 1
%             if abs(Normal_stress(i,j)) > sigma_t
%                 crack_path_1(i,j) = Normal_stress(i,j);
%             end
%         end
%     end
% end

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

node_configuration = Node_Configurator(Coordinates_Plate);
node_configuration = node_configuration';
node_configuration = flipud(node_configuration);
node_configuration = fliplr(node_configuration);

for k = 1:Number_of_Nodes %going over the nodes
    if crack_path_2(k,2) ~= 0 %chekcing if the point is overstressed
        Crack_Path_Draw = [];
        Length = 0;
        [row_node_configuration,column_node_configuration] = size(node_configuration);
        [row,column] = find(node_configuration == k);
        %Left Right Initiation
        left = column - 1;
        right = column + 1;
        no_more_right = false;
        no_more_left = false;

        %Up Down Initiation
        above = row - 1;
        below = row + 1;
        no_more_above = false;
        no_more_below = false;


        Directions = {};
        Orientations = {};

        Directions(1) = {[below,left]};
        Orientations(1) = {'W'};

        Directions(2) = {[row,left]};
        Orientations(2) = {'H'};

        Directions(3) = {[above,left]};
        Orientations(3) = {'D'};

        Directions(4) = {[below,column]};
        Orientations(4) = {'V'};

        Directions(5) = {[above,column]};
        Orientations(5) = {'V'};

        Directions(6) = {[below,right]};
        Orientations(6) = {'D'};

        Directions(7) = {[row,right]};
        Orientations(7) = {'H'};

        Directions(8) = {[above,right]};
        Orientations(8) = {'W'};

        for l = 1:8
            if Directions{l}(1)> 0 && Directions{l}(1) < row_node_configuration + 1 &&  Directions{l}(2)> 0 && Directions{l}(2) < column_node_configuration + 1

                if Orientations{l} ~= 'W'
                    Coordinates_Original = [Coordinates_Plate(node_configuration(row,column),1),Coordinates_Plate(node_configuration(row,column),2)];

                    Node_at_check = node_configuration(Directions{l}(1),Directions{l}(2));
                    Coordinates_Next = [Coordinates_Plate(Node_at_check,1),Coordinates_Plate(Node_at_check,2)];

                    X1 = Coordinates_Original(1);
                    X2 = Coordinates_Next(1);

                    Y1 =  Coordinates_Original(2);
                    Y2 =  Coordinates_Next(2);

                    Length = Length + (0.5) * sqrt(((X2 - X1)^2 + (Y2-Y1)^2));

                    Crack_Path_Temp = [Node_at_check, Coordinates_Next(1), Coordinates_Next(2), Orientations(l)];
                    Crack_Path_Draw = [Crack_Path_Draw; Crack_Path_Temp];
                end
            end
        end
        crack_path_2(k,3) = Length; %%assisning additions to the matrix
    end
end



end