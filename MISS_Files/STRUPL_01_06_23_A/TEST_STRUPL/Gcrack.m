function crack_path_2 = Gcrack(Normal_stress, nnd, connec, sigma_t, geom)

[rows,columns] = size(Normal_stress);
crack_path_1 = zeros(24,3);

Repitition_Remover = Node_Repitition_Remover(connec);

for i = 1:rows
    for j = 1:columns
        if Repitition_Remover(i,j) == 1
            if Normal_stress(i,j) > sigma_t
                crack_path_1(i,j) = Normal_stress(i,j);
            end
        end
    end
end


crack_path_2 = zeros(21,3);

%% Getting Element Numbers too

for k = 1:nnd
    for i = 1:rows
        for j = 1:columns
            if connec(i,j) == k 
                if Repitition_Remover(i,j) == 1
                    crack_path_2(k,1) = connec(i,j);
                    crack_path_2(k,2) = crack_path_1(i,j);
                end
            end
        end
    end
end


for k = 1:nnd %going over the nodes
    if crack_path_2(k,2) ~= 0 %chekcing if the point is overstressed
        Length = 0;
        for l = k-3:k+3 %Looping over nearby nodes
            if l > 1 && l < nnd
            if crack_path_2(l,2) ~= 0 %Checking if the nearby nodes are overstressed
                    % Finding Differences and adding them up
                    Length = Length + (0.5) * sqrt((geom(k,1) - geom(l,1))^2 +...
                            (geom(k,2) - geom(l,2))^2); 
            end
            end
        end
        crack_path_2(k,3) = Length; %%assisning additions to the matrix
    end
end



end