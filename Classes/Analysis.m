classdef Analysis

    %% Properties and Value Types
    properties
        number_of_dof_per_node = 0;
        dim = 0;
        Elastic_Modulus = 0;
        poissons_ratio = 0;
        fm = 0; % compressive strenght of mortar
        ft =0; % tensile strength of mortar

        % Fracture Parameters
        sigma_t = 0;
        Percentage_Limit_Tension = 0;
        tau = 0;
        Percentage_Limit_Shear = 0;
        Percentage_Fracture_Energy_Traction = 0;
        d_u = 0;
        Percentage_Ultimate_Deformation_Maximum_Step = 0;
        tol = 0;

        % Origin
        X_origin = 0;
        Y_origin = 0;

        % Gauss Points
        ngpb = 0;
        ngps = 0;

        % Structure Object
        STRUCTURE;

        is_dx_enabled = true;
        is_dy_enabled = true;
        is_dz_enabled = true;
        is_rx_enabled = true;
        is_ry_enabled = true;
        is_rz_enabled = true;

    end

    %% Methods
    methods
        function element = Element(Length_of_Element,Width_of_Element,number_of_degree_of_freedom_per_element,total_number_of_elements_connected,number_of_nodes_in_the_element)
            element.Length_of_Element = Length_of_Element;
            element.Width_of_Element = Width_of_Element;
            element.number_of_degree_of_freedom_per_element = number_of_degree_of_freedom_per_element;
            element.total_number_of_elements_connected = total_number_of_elements_connected;
            element.number_of_nodes_in_the_element = number_of_nodes_in_the_element;
        end
    end
end