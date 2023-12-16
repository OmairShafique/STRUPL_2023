classdef Node

    %% Properties and Value Types
    properties
        id = 0; % node identification number
        coord = []; % coordinates of node
        boundary_conditions = []; % essential boundary condition flags vector
        nodalLoad = []; % applied load components vector [fx fy fz mx my mz]
        external_load = []; % externally applied load on the node in matrix form
    end

    %% Methods
    methods
        % Constuctor Method
        function node = Node(id,coord,nodalLoad)
            node.id = id;
            node.coord = coord;
            node.nodalLoad = nodalLoad;
        end
    end


end










