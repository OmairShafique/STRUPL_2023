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
hold

%% Utilizing the Crack_path_matrix to draw the crack paths

Crack_Path_Matrix_Max = sortrows(Crack_Path,2);

node_configuration = Node_Configurator(nodal_coordinte_values);


[rows,columns] = size(Crack_Path_Matrix_Max);
% node_configuration = node_configuration';
% node_configuration = fliplr(node_configuration);
% node_configuration = flipud(node_configuration);

[rows_nodes_configurtion,columns_node_configuration] = size(node_configuration);


%% Using the highest Cracks from Crack_Path_matrix_max and making paths to the
% nearest horizontal vertical etc
for i = 1:rows
    % if Crack_Path_Matrix_Max(i,2) > 0
    %     break;
    % end

    [row,column] = find(node_configuration == Crack_Path_Matrix_Max(i,1));
    Crack_Path_Draw = [Crack_Path_Matrix_Max(i,1)];

    switch Direction
        case 'V' %Code for Horizontal
                left = column - 1;
                right = column + 1;
                no_more_right = false;
                no_more_left = false;
            while true
                if left > 0 && right < columns_node_configuration + 1
                    [row_for_left,c] = find(Crack_Path_Matrix_Max == node_configuration(row,left));
                    [row_for_right,c] = find(Crack_Path_Matrix_Max == node_configuration(row,right));

                    stress_value_left = Crack_Path_Matrix_Max(row_for_left,2);
                    stress_value_right = Crack_Path_Matrix_Max(row_for_right,2);

                    if stress_value_left > stress_value_right
                        if no_more_left
                            break;
                        end
                        % Now to Add Nodes to Drraw 
                        Crack_Path_Draw = [Crack_Path_Draw,node_configuration(row,left)];
                        left = left - 1;
                        no_more_right = true;                       
                    else
                        if no_more_right
                            break;
                        end
                        % Now to Add Nodes to Drraw
                        Crack_Path_Draw = [Crack_Path_Draw,node_configuration(row,right)];
                        right = right + 1;
                        no_more_left = true;                       
                    end
                else
                    break;
                end
            end
            
            Crack_Path_Draw = Crack_Path_Draw';
            
            [rows_crack,columns_crack] = size(Crack_Path_Draw);

            for k = 1:rows_crack
                Crack_Path_Draw(k,2) = nodal_coordinte_values(Crack_Path_Draw(k,1),1);
                Crack_Path_Draw(k,3) = nodal_coordinte_values(Crack_Path_Draw(k,1),2);
            end

            %Now to draw it
            
            plot(Crack_Path_Draw(:,2),Crack_Path_Draw(:,3),'Color','r');
            

        case 'H' %Code for Vertical
                above = row - 1;
                below = row + 1;
                no_more_above = false;
                no_more_below = false;
            while true
                if above > 0 && below < rows_nodes_configurtion + 1
                    [row_for_above,c] = find(Crack_Path_Matrix_Max == node_configuration(above,column));
                    [row_for_below,c] = find(Crack_Path_Matrix_Max == node_configuration(below,column));

                    stress_value_above = Crack_Path_Matrix_Max(row_for_above,2);
                    stress_value_below = Crack_Path_Matrix_Max(row_for_below,2);

                    if stress_value_above > stress_value_below
                        if no_more_above
                            break;
                        end
                        % Now to Add Nodes to Drraw 
                        Crack_Path_Draw = [Crack_Path_Draw,node_configuration(above,column)];
                        above = above - 1;
                        no_more_below = true;                       
                    else
                        if no_more_below
                            break;
                        end
                        % Now to Add Nodes to Drraw
                        Crack_Path_Draw = [Crack_Path_Draw,node_configuration(below,column)];
                        below = below + 1;
                        no_more_above = true;                       
                    end
                else
                    break;
                end
            end
            
            Crack_Path_Draw = Crack_Path_Draw';
            
            [rows_crack,columns_crack] = size(Crack_Path_Draw);

            for k = 1:rows_crack
                Crack_Path_Draw(k,2) = nodal_coordinte_values(Crack_Path_Draw(k,1),1);
                Crack_Path_Draw(k,3) = nodal_coordinte_values(Crack_Path_Draw(k,1),2);
            end
            %Now to draw it
            plot(Crack_Path_Draw(:,2),Crack_Path_Draw(:,3),'Color','b');
        
        case 'D' %Code for Diagonal


    end
end







end