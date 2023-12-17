classdef Analysis

    %% Properties and Value Types
    properties
        thickness_of_plate = 0; %thickness of the structure
        number_of_dof_per_node = 0;
        dim = 0;
        Elastic_Modulus = 0;
        poissons_ratio = 0;
        fm = 0; % compressive strenght of mortar
        ft =0; % tensile strength of mortar
        Boundary_Conditions = []; %boundary conditions matrix input directly in the form of a matrix
        External_Load = []; %external load applied on the structure specifying node and the dof associated
        gamma = 0; %weight gamma
        geom = []; %geometry of the strucutre matrix
        connec = []; %.
        nodal_coordinate_values = []; %coordinates of all the nodes
        nodal_connectivity_values = []; %element and node connectivity values




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
        function analysis = Analysis(thickness_of_plate,number_of_dof_per_node,dim,Boundary_Conditions,External_Load,nodal_coordinate_values,nodal_connectivity_values)
            analysis.thickness_of_plate = thickness_of_plate;
            analysis.number_of_dof_per_node = number_of_dof_per_node;
            analysis.dim = dim;
            analysis.Boundary_Conditions = Boundary_Conditions;
            analysis.External_Load = External_Load;

            analysis.STRUCTURE = Structure(nodal_coordinate_values,nodal_connectivity_values);
        end


    end
end