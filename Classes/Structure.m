    classdef Structure
        properties
            type_of_element = []
            thickness_of_plate = 0; %thickness of the structure
            Boundary_Conditions = [];
            External_Load = [];
            gamma = 0;
            nodal_coordinate_values = [];
            nodal_connectivity_values = [];
            geom = [];
            connec = [];

            Element_Type = [];

            Number_of_Elements = length(nodal_connectivity_values);
            Number_of_Nodes = length(nodal_coordinate_values);
            Total_System_Degrees_of_Freedom = Number_of_Nodes * number_of_dof_per_node;
            Degrees_of_Freedom_Per_Element = number_of_nodes_per_element * number_of_dof_per_node;

        end
    end