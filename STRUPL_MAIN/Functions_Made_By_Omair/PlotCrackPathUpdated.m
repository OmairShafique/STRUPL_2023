function PlotCrackPathUpdated(Cracks_All,Crack_Path,nodal_coordinte_values,nodal_connectivity_values)

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

plot(X,Y,'k')
fill(X,Y,'w')
set(f1,'name','Crack Path','numbertitle','off');
title('Crack Path Estimate');
axis off ;
hold

%% Utilizing the Crack_path_matrix to draw the crack paths

% Crack_Path_Matrix_Max = sortrows(Crack_Path,2);

node_configuration = Node_Configurator(nodal_coordinte_values);

node_configuration = node_configuration';
node_configuration = flipud(node_configuration);
node_configuration = fliplr(node_configuration);

[rows,columns] = size(Crack_Path);


[rows_nodes_configurtion,columns_node_configuration] = size(node_configuration);

%% Using the highest Cracks from Crack_Path_matrix_max and making paths to the
% nearest horizontal vertical etc
for i = 1:rows
    if Crack_Path(i) ~= 0
        [rows_e,columns_e] = size(Cracks_All{i});
        [row,column] = find(node_configuration == Crack_Path(i,1));
        Coordinates_Original = [nodal_coordinte_values(node_configuration(row,column),1),nodal_coordinte_values(node_configuration(row,column),2)];

        for j = 1:rows_e
            Coordinates_New = [Cracks_All{i}{j,2},Cracks_All{i}{j,3}];

            X = [Coordinates_Original(1,1), Coordinates_New(1,1)];
            Y = [Coordinates_Original(1,2), Coordinates_New(1,2)];

            Orientation = Cracks_All{i}{j,4};

            switch Orientation
                case 'H'
                    plot(X,Y,'Color','r');
                case 'V'
                    plot(X,Y,'Color','b');
                case 'D'
                    plot(X,Y,'Color','g');
                case 'W'
                    % plot(X,Y,'Color','black');
            end
        end


    end
end
end