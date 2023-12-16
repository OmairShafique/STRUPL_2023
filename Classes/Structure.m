classdef Structure

    %% Properties and Value Types
    properties
        %% Input Values
        thickness_of_plate = 0; %thickness of the structure
        Boundary_Conditions = []; %boundary conditions matrix input directly in the form of a matrix
        External_Load = []; %external load applied on the structure specifying node and the dof associated
        gamma = 0; %weight gamma
        nodal_coordinate_values = []; %coordinates of all the nodes
        nodal_connectivity_values = []; %element and node connectivity values
        geom = []; %geometry of the strucutre matrix
        connec = []; %

        Element_Type = []; % Element Type i.e. Triangular, Square etc

        ELEMENTS;

        %% Infered Values
        Number_of_Elements = length(nodal_connectivity_values); % Number of Elements in the Structure
        Number_of_Nodes = length(nodal_coordinate_values); % Number of nodes in the structure
        Total_System_Degrees_of_Freedom = Number_of_Nodes * number_of_dof_per_node; % Total System Degrees of Freedom calculated
        Degrees_of_Freedom_Per_Element = number_of_nodes_per_element * number_of_dof_per_node; % Calculating the Degrees of Freedom Per Element

    end

    %% Methods
    methods
        % Constuctor Method
        function structure = Structure(thickness_of_plate,Boundary_Conditions,External_Load,gamma,nodal_coordinate_values,nodal_connectivity_values,geom,connec)
            structure.thickness_of_plate = thickness_of_plate;
            structure.Boundary_Conditions = Boundary_Conditions;
            structure.External_Load = External_Load;
            structure.gamma = gamma;
            structure.nodal_coordinate_values = nodal_coordinate_values;
            structure.nodal_connectivity_values = nodal_connectivity_values;
            structure.geom = geom;
            structure.connec = connec;
        end
    end
end