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
        % Constructor
        function analysis = Analysis(thickness_of_plate,number_of_dof_per_node,dim,Boundary_Conditions,External_Load,nodal_coordinate_values,nodal_connectivity_values)
            analysis.thickness_of_plate = thickness_of_plate;
            analysis.number_of_dof_per_node = number_of_dof_per_node;
            analysis.dim = dim;
            analysis.Boundary_Conditions = Boundary_Conditions;
            analysis.External_Load = External_Load;

            analysis.STRUCTURE = Structure(nodal_coordinate_values,nodal_connectivity_values);
        end

        % Populating nf
        function [nf,total_numbers_of_active_dof,nf_g] = Populating_nf()
            %%TEMPORARY SPACE Populating the NF matrix-------------
            nf = ones(Number_of_Nodes,this.number_of_dof_per_node); %nodal freedom matrix set to zeros

            node_number_where_active = this.Boundary_Conditions(:,1); %creating a single coumn matrix with node numbers where dof are released
            dof_x_displacement = this.Boundary_Conditions(:,2); %dof release matrix from x-disp column in Boundry condition file
            dof_y_displacement = this.Boundary_Conditions(:,3); %dof release matrix from y-disp column in Boundry condition file
            dof_z_displacement = this.Boundary_Conditions(:,4); %dof release matrix from z-disp column in Boundry condition file
            dof_x_rotation = this.Boundary_Conditions(:,5); %dof release matrix from x-rotation column in Boundry condition file
            dof_y_rotation = this.Boundary_Conditions(:,6); %dof release matrix from y-rotation column in Boundry condition file
            dof_z_rotation = this.Boundary_Conditions(:,7); %dof release matrix from z-rotation column in Boundry condition file
            %miss = Boundry_Conditions(:,8); %dof release matrix from last column in Boundry condition file

            number_of_active_boundry_conditions = length(node_number_where_active); %total number of nodes at which releases present


            %meant to populate the nf matrix
            for i = 1:number_of_active_boundry_conditions
                if(this.dim == 2) %checking if 2D
                    if this.is_dx_enabled
                        nf(node_number_where_active(i),1) = dof_x_displacement(i); %Changing nf at x_displacmeent
                    end
                    if this.is_dy_enabled
                        nf(node_number_where_active(i),2) = dof_y_displacement(i);
                    end
                    if this.is_rx_enabled
                        nf(node_number_where_active(i),3) = dof_x_rotation(i);
                    end
                else
                    if(this.dim ==3) %Checking if 3D
                        if this.is_dx_enabled
                            nf(node_number_where_active(i),1) = dof_x_displacement(i); %Changing nf at x_displacmeent
                        end
                        if this.is_dy_enabled
                            nf(node_number_where_active(i),2) = dof_y_displacement(i); %Changing nf at x_displacmeent
                        end
                        if this.is_dz_enabled
                            nf(node_number_where_active(i),3) = dof_z_displacement(i); %Changing nf at x_displacmeent
                        end
                        if this.is_rx_enabled
                            nf(node_number_where_active(i),4) = dof_x_rotation(i); %Changing nf at x_displacmeent
                        end
                        if this.is_ry_enabled
                            nf(node_number_where_active(i),5) = dof_y_rotation(i); %Changing nf at x_displacmeent
                        end
                        if this.is_rz_enabled
                            nf(node_number_where_active(i),6) = dof_z_rotation(i); %Changing nf at x_displacmeent
                        end
                    end
                end
            end

            total_numbers_of_active_dof = 0;

            for i = 1:length(nf)
                for j = 1:this.number_of_dof_per_node
                    if(nf(i,j) == 1) %checking if nf has 1 or 0
                        total_numbers_of_active_dof = total_numbers_of_active_dof + 1; %if 1 then increase value by 1
                    end
                end
            end


            nf_g = nf;
            n = 0;
            for i=1:Number_of_Nodes
                for j=1:this.number_of_dof_per_node
                    if nf_g(i,j) ~= 0
                        n=n+1;
                        nf_g(i,j)=n;
                    end
                end
            end

        end

        %Populating Load Matric
        function [Load,Global_force_vector,Gravity_Load] = Populating_Load(nf_g)
            Nodal_load = this.External_load(:,1);
            Force = this.External_load(:,2:4);

            % 1.0 Assign Concentrated load  is the first step and it is ok%
            Load = zeros(Number_of_Nodes,3);

            for i=1:length(Nodal_load)
                Load(Nodal_load(i),1:3) = Force(i,1:3);
            end

            Global_force_vector = zeros(total_numbers_of_active_dof,1);
            for i=1:Number_of_Nodes
                for j=1:this.number_of_dof_per_node
                    if nf_g(i,j) ~= 0
                        Global_force_vector(nf_g(i,j))= Load(i,j);
                    end
                end
            end

            % 2.0 Assign gravity load that generates body forces
            Gravity_Load=[0 ,-this.gamma]';
        end

        %



    end
end