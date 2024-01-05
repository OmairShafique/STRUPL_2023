classdef Structure

    %% Properties and Value Types
    properties
        %% Input Values
        Length = 0; %Assumed constant length of element
        Width = 0; %Assumed constatn width of element
        nodal_coordinate_values = []; %coordinates of all the nodes
        nodal_connectivity_values = []; %element and node connectivity values

        number_of_nodes_per_element = 0;

        ELEMENTS;

        %% Infered Values
        Number_of_Elements = 0; % Number of Elements in the Structure
        Number_of_Nodes = 0; % Number of nodes in the structure
        Total_System_Degrees_of_Freedom = 0; % Total System Degrees of Freedom calculated
        Degrees_of_Freedom_Per_Element = 0; % Calculating the Degrees of Freedom Per Element
    end

    %% Methods
    methods
        % Constuctor Method
        function structure = Structure(Length,Width,nodal_coordinate_values,nodal_connectivity_values,number_of_nodes_per_element)
            structure.Length = Length;
            structure.Width = Width;
            structure.nodal_coordinate_values = nodal_coordinate_values;
            structure.nodal_connectivity_values = nodal_connectivity_values;
            structure.number_of_nodes_per_element = number_of_nodes_per_element;

            structure.Number_of_Elements = length(nodal_connectivity_values); % Number of Elements in the Structure
            structure.Number_of_Nodes = length(nodal_coordinate_values); % Number of nodes in the structure
            % structure.Total_System_Degrees_of_Freedom = structure.Number_of_Nodes * structure.number_of_dof_per_node; % Total System Degrees of Freedom calculated
            % structure.Degrees_of_Freedom_Per_Element = number_of_nodes_per_element * number_of_dof_per_node; % Calculating the Degrees of Freedom Per Element

            Elements = Element(1,structure.Number_of_Elements,Node(1,1));
            for i = 1:structure.Number_of_Elements
                nodes = Node(1,number_of_nodes_per_element);
                for j = 1:number_of_nodes_per_element
                    coord = structure.nodal_coordinate_values(structure.nodal_connectivity_values(i,j),:);
                    id = structure.nodal_connectivity_values(i,j);
                    node = Node(id,coord);
                    nodes(j) = node;
                end
                element = Element(structure.Length,structure.Width,nodes);
                Elements(i) = element;
            end
            structure.ELEMENTS = Elements';
            
        end
    end
end