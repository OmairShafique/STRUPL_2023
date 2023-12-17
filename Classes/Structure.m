classdef Structure

    %% Properties and Value Types
    properties
        %% Input Values
        Length = 0; %Assumed constant length of element
        Width = 0; %Assumed constatn width of element
        nodal_coordinate_values = []; %coordinates of all the nodes
        nodal_connectivity_values = []; %element and node connectivity values

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
        function structure = Structure(nodal_coordinate_values,nodal_connectivity_values)
            structure.nodal_coordinate_values = nodal_coordinate_values;
            structure.nodal_connectivity_values = nodal_connectivity_values;
        end

        function this.PopulatingElements()
            Elements = [1,this.Number_of_Elements];
            for i = 1:this.Number_of_Elements
                nodes = [1,number_of_nodes_per_element];
                    for j = 1:number_of_nodes_per_element
                        coord = this.nodal_coordinate_values(this.nodal_connectivity_values(i,j));
                        id = this.nodal_connectivity_values(i,j);

                        node = Node(id,coord);
                        nodes(j) = node;
                    end
                    element = Element(this.Length,this.Width,nodes);
                    Elements(i) = element;
            end
            this.ELEMENTS = Elements;
        end
    end
end