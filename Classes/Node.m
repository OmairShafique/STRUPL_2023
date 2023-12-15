classdef Node

    properties
        id = 0; % node identification number 
        coord = []; % coordinates of node
        ebc = []; % essential boundary condition flags vector
        nodalLoad = []; % applied load components vector [fx fy fz mx my mz]  
        external_load = 0;
        prescDispl = []; % prescribed displacement values vector [dx dy dz rx ry rz]
       


    end

    methods
        %Constuctor Method
        function node = Node(id,coord,ebc,nodalLoad,prescDispl)
            node.id = id;
            node.coord = coord;
            node.ebc = ebc;
            node.nodalLoad = nodalLoad;
            node.prescDispl= prescDispl;

        end


    end


end










