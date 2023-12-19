function nodes_configuration = Node_Configurator(Coordinates)
%% Going over the coordinates file and creating a matrix with the node network so it can be used later to find the crack path


[rows,columns] = size(Coordinates);
nodes_configuration = 0;

% %% Checking to see which column is changing
% Column_Changing = 1;
% 
% if Coordinates(1,1) ~= Coordinates(2,1) %if the first column is being varied, 
%     Column_Changing = 2; %It means the second column is the pivot
% else
%     Column_Changing = 1; %Or otherwise
% end


Row_Check = 0;
current_Row = 1;
row_index = 2;

Xs = Coordinates(:,1);
Ys = Coordinates(:,2);

max_X = max(Xs);
max_Y = max(Ys);

min_X = min(Xs);
min_Y = min(Ys);


Range_X = max_X - min_X;
Range_Y = max_Y - min_Y;

X_sorted = sort(Xs);
Y_sorted = sort(Ys);


numb_of_rows = 1;
numb_of_columns = 1;

X_Values = [min_X];
Y_Values = [min_Y];

for i = 2:size(X_sorted)
    if X_sorted(i) ~= X_sorted(i-1)
        numb_of_rows = numb_of_rows + 1;
        X_Values = [X_Values;X_sorted(i)];
    end
end

for i = 2:size(Y_sorted)
    if Y_sorted(i) ~= Y_sorted(i-1)
        numb_of_columns = numb_of_columns + 1;
        Y_Values = [Y_Values;Y_sorted(i)];
    end
end

    new_coordinates = zeros(rows,columns);

    for i = 1:rows
            row = find(X_Values == Coordinates(i,1));
            column = find(Y_Values == Coordinates(i,2));
            new_coordinates(i,1) = row;
            new_coordinates(i,2) = column;
    end
    
    nodes_configuration = zeros(numb_of_rows,numb_of_columns);

for i = 1:rows
    nodes_configuration(new_coordinates(i,1),new_coordinates(i,2)) = i;
end


end