function Cracks_All =  CrackPropagation(Crack_Path,nodal_coordinte_values,nodal_connectivity_values)

%% Drawings the Graph and the Elements with Nodes
nel = length(nodal_connectivity_values);
nnd = length(nodal_coordinte_values);
nne = size(nodal_connectivity_values,2);

X = zeros(nne,nel);
Y = zeros(nne,nel);

for iel=1:nel
    for i=1:nne
        nd(i)=nodal_connectivity_values(iel,i);      % extract connected node for (iel)-th element
        X(i,iel)=nodal_coordinte_values(nd(i),1);    % extract x value of the node
        Y(i,iel)=nodal_coordinte_values(nd(i),2);    % extract y value of the node
    end
end

% Plotting the FEM mesh, diaplay Node numbers and Element numbers


f1 = figure ;
set(f1,'name','Crack Path','numbertitle','off') ;
plot(X,Y,'k')
fill(X,Y,'w')
title('Crack Path Estimate') ;



axis off ;
hold

%% Utilizing the Crack_path_matrix to draw the crack paths

Cracks_All = {};
Crack_Path_Matrix_Max = sortrows(Crack_Path,2);


node_configuration = Node_Configurator(nodal_coordinte_values);



[rows,columns] = size(Crack_Path_Matrix_Max);


[rows_nodes_configurtion,columns_node_configuration] = size(node_configuration);



%% Using the highest Cracks from Crack_Path_matrix_max and making paths to the
% nearest horizontal vertical etc
for i = 1:rows

    % if Crack_Path_Matrix_Max(i,2) > 0
    %     break;
    % end

    [row,column] = find(node_configuration == Crack_Path_Matrix_Max(i,1));
    Crack_Path_Draw = [];
    Higest_Node_Stress = Crack_Path_Matrix_Max(i,1);

    while true

        %Left Right Initiation
        left = column - 1;
        right = column + 1;
        % no_more_right = false;
        % no_more_left = false;

        %Up Down Initiation
        above = row - 1;
        below = row + 1;
        % no_more_above = false;
        % no_more_below = false;

        if left < 1 || right > columns_node_configuration+1 || above < 1 || below > rows_nodes_configurtion+1
            break;
        end


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

        %If you want to change the code criteria you have to change this
        %part
        This_Node = Higest_Node_Stress;
        for l = 1:8
            if Directions{l}(1)> 0 && Directions{l}(1) < rows_nodes_configurtion + 1 &&  Directions{l}(2)> 0 && Directions{l}(2) < columns_node_configuration + 1
                Node_at_check = node_configuration(Directions{l}(1),Directions{l}(2));
                if Node_at_check ~= 0
                    if Crack_Path(Node_at_check,2) ~= 0
                        if abs(Crack_Path(Higest_Node_Stress,2)) < abs(Crack_Path(Node_at_check,2))
                            Higest_Node_Stress = Node_at_check;
                        end
                    end
                end

            end
        end

        if Higest_Node_Stress == This_Node
            break;
        end

        Coordinates_Original = [nodal_coordinte_values(This_Node,1),nodal_coordinte_values(This_Node,2)];
        Coordinates_Next = [nodal_coordinte_values(Higest_Node_Stress,1),nodal_coordinte_values(Higest_Node_Stress,2)];
        Crack_Path_Temp = [Higest_Node_Stress, Coordinates_Next(1), Coordinates_Next(2)];
        Crack_Path_Draw = [Crack_Path_Draw; Crack_Path_Temp];

        [row,column] = find(node_configuration == Higest_Node_Stress);
    end

    if size(Crack_Path_Draw) > 1
        plot(Crack_Path_Draw(:,2),Crack_Path_Draw(:,3),'Color','r');
        [row_Cracks,column_Cracks] = size(Cracks_All);
        Cracks_All(row_Cracks + 1,1) = {Crack_Path_Draw};
    end


end


end