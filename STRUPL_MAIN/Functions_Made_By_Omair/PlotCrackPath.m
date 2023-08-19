function PlotCrackPath(Crack_Path,nodal_coordinte_values,nodal_connectivity_values,Direction)

%% Drawings the Graph and the Elements with Nodes
nel = length(nodal_connectivity_values);
nnd = length(nodal_coordinte_values);
nne = size(nodal_connectivity_values,2);

X = zeros(nne,nel);
Y = zeros(nne,nel);

for iel=1:nel
    for i=1:nne
        nd(i)=nodal_connectivity_values(iel,i);         % extract connected node for (iel)-th element
        X(i,iel)=nodal_coordinte_values(nd(i),1);    % extract x value of the node
        Y(i,iel)=nodal_coordinte_values(nd(i),2);    % extract y value of the node
    end
end

% Plotting the FEM mesh, diaplay Node numbers and Element numbers
f1 = figure ;
set(f1,'name','Mesh','numbertitle','off') ;
plot(X,Y,'k')
fill(X,Y,'w')

title('Finite Element Mesh') ;
axis off ;

%% Utilizing the Crack_path_matrix to draw the crack paths

Crack_Path_Matrix_Max = sortrows(Crack_Path,2);

node_configuration = Node_Configurator(nodal_coordinte_values);


[rows,columns] = size(Crack_Path_Matrix_Max);
[rows_nodes_configurtion,columns_node_configuration] = size(node_configuration);


%% Using the highest Cracks from Crack_Path_matrix_max and making paths to the
% nearest horizontal vertical etc
for i = 1:rows
    [row,column] = find(node_configuration == Crack_Path_Matrix_Max(i,1));

    switch Direction
        case 'Horizontal'
            while true
                left = row - 1;
                right = row + 1;
                no_more_right = false;
                no_more_left = false;

                if left < 1 || right > rows_nodes_configurtion
                    [row_for_left,c] = find(Crack_Path_Matrix_Max == node_configuration(left,column));
                    [row_for_right,c] = find(Crack_Path_Matrix_Max == node_configuration(right,column));

                    stress_value_left = Crack_Path_Matrix_Max(row_for_left,2);
                    stress_value_right = Crack_Path_Matrix_Max(row_for_right,2);

                    if stress_value_left > stress_value_right
                        left = left - 1;
                        no_more_right = true;
                    else
                        right = right + 1;
                        no_more_left = true;
                    end

                else
                    break;
                end
            end

        case 'Vertical'

        case 'Diagonal'

    end
end







end