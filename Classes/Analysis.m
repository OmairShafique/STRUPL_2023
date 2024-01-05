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
        geom = []; %analysisObject.geometry of the strucutre matrix
        connec = []; %.
        ngps = 0; % Number of Gauss Points Sheer
        ngpb = 0; % Number of Gauss Points Bending
        Length_of_Element = 0;
        Width_of_Element = 0;
        Element_Type = 0; % Element Type i.e. Triangular, Square etc
        nf_g = 0;


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
        function analysis = Analysis(Length_of_Element,Width_of_Element ...
                ,thickness_of_plate,number_of_dof_per_node,dim,Boundary_Conditions ...
                ,External_Load,nodal_coordinate_values,nodal_connectivity_values,Element_Type)
            analysis.thickness_of_plate = thickness_of_plate;
            analysis.number_of_dof_per_node = number_of_dof_per_node;
            analysis.dim = dim;
            analysis.Boundary_Conditions = Boundary_Conditions;
            analysis.External_Load = External_Load;
            analysis.Length_of_Element = Length_of_Element;
            analysis.Width_of_Element = Width_of_Element;
            analysis.geom = nodal_coordinate_values;
            analysis.connec = nodal_connectivity_values;
            analysis.Element_Type = Element_Type;

            analysis.STRUCTURE = Structure(Length_of_Element,Width_of_Element,nodal_coordinate_values,nodal_connectivity_values,Element_Type); 
            analysis.STRUCTURE.Degrees_of_Freedom_Per_Element = Element_Type * number_of_dof_per_node;
        end
    end

    methods (Static)
        % Main Class where all logic will follow
        function Engine(analysisObject)
            %% POPULATING NF
            nf = ones(analysisObject.STRUCTURE.Number_of_Nodes,analysisObject.number_of_dof_per_node); %nodal freedom matrix set to zeros

            node_number_where_active = analysisObject.Boundary_Conditions(:,1); %creating a single coumn matrix with node numbers where dof are released
            dof_x_displacement = analysisObject.Boundary_Conditions(:,2); %dof release matrix from x-disp column in Boundry condition file
            dof_y_displacement = analysisObject.Boundary_Conditions(:,3); %dof release matrix from y-disp column in Boundry condition file
            dof_z_displacement = analysisObject.Boundary_Conditions(:,4); %dof release matrix from z-disp column in Boundry condition file
            dof_x_rotation = analysisObject.Boundary_Conditions(:,5); %dof release matrix from x-rotation column in Boundry condition file
            dof_y_rotation = analysisObject.Boundary_Conditions(:,6); %dof release matrix from y-rotation column in Boundry condition file
            dof_z_rotation = analysisObject.Boundary_Conditions(:,7); %dof release matrix from z-rotation column in Boundry condition file
            %miss = Boundry_Conditions(:,8); %dof release matrix from last column in Boundry condition file

            number_of_active_boundry_conditions = length(node_number_where_active); %total number of nodes at which releases present


            %meant to populate the nf matrix
            for i = 1:number_of_active_boundry_conditions
                if(analysisObject.dim == 2) %checking if 2D
                    if analysisObject.is_dx_enabled
                        nf(node_number_where_active(i),1) = dof_x_displacement(i); %Changing nf at x_displacmeent
                    end
                    if analysisObject.is_dy_enabled
                        nf(node_number_where_active(i),2) = dof_y_displacement(i);
                    end
                    if analysisObject.is_rx_enabled
                        nf(node_number_where_active(i),3) = dof_x_rotation(i);
                    end
                else
                    if(analysisObject.dim ==3) %Checking if 3D
                        if analysisObject.is_dx_enabled
                            nf(node_number_where_active(i),1) = dof_x_displacement(i); %Changing nf at x_displacmeent
                        end
                        if analysisObject.is_dy_enabled
                            nf(node_number_where_active(i),2) = dof_y_displacement(i); %Changing nf at x_displacmeent
                        end
                        if analysisObject.is_dz_enabled
                            nf(node_number_where_active(i),3) = dof_z_displacement(i); %Changing nf at x_displacmeent
                        end
                        if analysisObject.is_rx_enabled
                            nf(node_number_where_active(i),4) = dof_x_rotation(i); %Changing nf at x_displacmeent
                        end
                        if analysisObject.is_ry_enabled
                            nf(node_number_where_active(i),5) = dof_y_rotation(i); %Changing nf at x_displacmeent
                        end
                        if analysisObject.is_rz_enabled
                            nf(node_number_where_active(i),6) = dof_z_rotation(i); %Changing nf at x_displacmeent
                        end
                    end
                end
            end

            total_numbers_of_active_dof = 0;

            for i = 1:length(nf)
                for j = 1:analysisObject.number_of_dof_per_node
                    if(nf(i,j) == 1) %checking if nf has 1 or 0
                        total_numbers_of_active_dof = total_numbers_of_active_dof + 1; %if 1 then increase value by 1
                    end
                end
            end


            analysisObject.nf_g = nf;
            n = 0;
            for i=1:analysisObject.STRUCTURE.Number_of_Nodes
                for j=1:analysisObject.number_of_dof_per_node
                    if analysisObject.nf_g(i,j) ~= 0
                        n=n+1;
                        analysisObject.nf_g(i,j)=n;
                    end
                end
            end

            clear node_number_where_active
            clear dof_x_displacement
            clear dof_y_displacement
            clear dof_z_displacement
            clear dof_x_rotation
            clear dof_y_rotation
            clear dof_z_rotation
            clear miss

            %%
            %% LOADING
            Nodal_load = analysisObject.External_Load(:,1);
            Force = analysisObject.External_Load(:,2:4);

            % 1.0 Assign Concentrated load  is the first step and it is ok%
            Load = zeros(analysisObject.STRUCTURE.Number_of_Nodes,3);

            for i=1:length(Nodal_load)
                Load(Nodal_load(i),1:3) = Force(i,1:3);
            end

            Global_force_vector = zeros(total_numbers_of_active_dof,1);
            for i=1:analysisObject.STRUCTURE.Number_of_Nodes
                for j=1:analysisObject.number_of_dof_per_node
                    if analysisObject.nf_g(i,j) ~= 0
                        Global_force_vector(analysisObject.nf_g(i,j))= Load(i,j);
                    end
                end
            end

            % 2.0 Assign gravity load that generates body forces
            % fg_gravity = fg_matrix_calculator(analysisObject);

            deeb=formdeeb(analysisObject.Elastic_Modulus,analysisObject.poissons_ratio,analysisObject.thickness_of_plate); % Matrix of elastic properties for plate bending
            dees=formdees(analysisObject.Elastic_Modulus,analysisObject.poissons_ratio,analysisObject.thickness_of_plate); % Matrix of elastic properties for plate shear
            %%
            %% STIFFNESS MATRIX
            % 1.Form the matrix containing the abscissas and the weights of Gauss points
            sampb = gauss(analysisObject.ngpb);
            samps = gauss(analysisObject.ngps);

            % 2. initialize the global stiffness matrix to zero
            Global_stiffness_matrix = zeros(total_numbers_of_active_dof,total_numbers_of_active_dof);

            for iel = 1:analysisObject.STRUCTURE.Number_of_Elements
                if analysisObject.Element_Type == 3 && analysisObject.ngpb == 0
                    [bee,fun_3,g,A,d_3] = elem_T3(iel,analysisObject);
                    dee=formdeeb(analysisObject.Elastic_Modulus,analysisObject.poissons_ratio,analysisObject.thickness_of_plate);
                    ke=analysisObject.thickness_of_plate*A*bee'*dee*bee; % Integrate stiffness matrix
                    Global_stiffness_matrix=form_KK(Global_stiffness_matrix,ke, g);
                end
                [coord, g] = platelem_q4(iel,analysisObject); % coordinates of the nodes of element i,
                % and its steering vector
                keb = zeros(analysisObject.STRUCTURE.Degrees_of_Freedom_Per_Element,analysisObject.STRUCTURE.Degrees_of_Freedom_Per_Element) ; % Initialize the element bending
                % stiffness matrix to zero
                kes=zeros(analysisObject.STRUCTURE.Degrees_of_Freedom_Per_Element,analysisObject.STRUCTURE.Degrees_of_Freedom_Per_Element) ; % Initialize the element Shear
                % stiffness matrix to zero

                %--------------------------------------------------------------------------
                % Integrate element Shear stiffness and assemble it in global matrix
                %--------------------------------------------------------------------------
                for ig=1: analysisObject.ngps
                    wi = samps(ig,2);
                    for jg=1: analysisObject.ngps
                        wj=samps(jg,2);
                        [der,fun] = fmquad(samps, analysisObject, ig,jg); % Derivative of shape functions
                        % in local coordinates
                        jac=der'*coord;                    % Compute Jacobian matrix
                        d=det(jac);                       % Compute determinant of
                        % Jacobian matrix
                        jac1=inv(jac);                    % Compute inverse of the
                        % Jacobian
                        deriv=jac1*der';                   % Derivative of shape functions
                        % in global coordinates
                        bees=formbees(deriv,fun,analysisObject.STRUCTURE.number_of_nodes_per_element,analysisObject.STRUCTURE.Degrees_of_Freedom_Per_Element,analysisObject.number_of_dof_per_node); % Form matrix [B]
                        kes=kes + (5/6)*d*wi*wj*bees'*dees*bees; % Integrate stiffness matrix
                    end
                end
                Global_stiffness_matrix=form_KK(Global_stiffness_matrix,kes, g,analysisObject.STRUCTURE.Degrees_of_Freedom_Per_Element); % assemble global stiffness matrix
            end

            %%
            %% DISPLAYING RESULTS

            delta = Global_stiffness_matrix\Global_force_vector;

            if analysisObject.Element_Type == 3 %Parlas
                for i=1: analysisObject.STRUCTURE.Number_of_Nodes %
                    if analysisObject.nf_g(i,1) == 0 %
                        x_disp =0.; %
                    else
                        x_disp = delta(analysisObject.nf_g(i,1)); %
                    end
                    %
                    if analysisObject.nf_g(i,2) == 0 %
                        y_disp = 0.; %
                    else
                        y_disp = delta(analysisObject.nf_g(i,2)); %
                    end
                    disp([i x_disp y_disp]) % Display displacements of each node
                    DISP(i,:) = [ x_disp y_disp];
                end
            end

            format short e
            disp('node w_disp x_slope y_slope ') %
            for i=1: analysisObject.STRUCTURE.Number_of_Nodes %
                if analysisObject.nf_g(i,1) == 0 %
                    w_disp =0.; %
                else
                    w_disp = delta(analysisObject.nf_g(i,1)); %
                end
                %
                if analysisObject.nf_g(i,2) == 0 %
                    x_slope = 0.; %
                else
                    x_slope = delta(analysisObject.nf_g(i,2)); %
                end
                %
                if analysisObject.number_of_dof_per_node > 2
                    if analysisObject.nf_g(i,3) == 0 %
                        y_slope = 0.; %
                        disp([i w_disp x_slope y_slope]) % Display displacements of each node
                        DISP(i,:) = [ w_disp x_slope y_slope];
                    else
                        y_slope = delta(analysisObject.nf_g(i,3)); %
                        disp([i w_disp x_slope y_slope]) % Display displacements of each node
                        DISP(i,:) = [ w_disp x_slope y_slope];
                    end
                else
                    disp([i w_disp x_slope]) % Display displacements of each node
                    DISP(i,:) = [ w_disp x_slope];
                end
            end

            disp('Calculate moments and shear forces the center of each element')
            ngp=1; %
            samp=gauss(ngp);
            %
            % Still needs Implementation
            for i=1:analysisObject.STRUCTURE.Number_of_Elements
                [coord,g] = platelem_q4(i,analysisObject); % coordinates of the nodes of element i,
                % and its steering vector
                eld=zeros(analysisObject.STRUCTURE.Degrees_of_Freedom_Per_Element,1); % Initialize element displacement to zero
                for m=1:analysisObject.STRUCTURE.Degrees_of_Freedom_Per_Element %
                    if g(m)==0 %
                        eld(m)=0.; %
                    else %
                        eld(m)=delta(g(m)); % Retrieve element displacement from the
                        % global displacement vector
                    end
                end
                %
                for ig=1: ngp
                    wi = samp(ig,2);
                    for jg=1: ngp
                        wj=samp(jg,2);
                        [der,fun] = fmquad(samp, analysisObject, ig,jg,i); % Derivative of shape functions
                        % in local coordinates
                        jac=der'*coord; % Compute Jacobian matrix
                        d=det(jac); % Compute the determinant of
                        % Jacobian matrix
                        jac1=inv(jac); % Compute inverse of the Jacobian
                        deriv=jac1*der'; % Derivative of shape functions
                        % in global coordinates
                        %
                        beeb=formbeeb(deriv,analysisObject.STRUCTURE.number_of_nodes_per_element,analysisObject.STRUCTURE.Degrees_of_Freedom_Per_Element); % Form matrix [B_b]
                        chi_b = beeb*eld ; % compute bending curvatures
                        Moment = deeb*chi_b ; % Compute moments
                        bees=formbees(deriv,fun,analysisObject.STRUCTURE.number_of_nodes_per_element,analysisObject.STRUCTURE.Degrees_of_Freedom_Per_Element,analysisObject.number_of_dof_per_node); % Form matrix [B_s]
                        chi_s = bees*eld ; % compute shear curvatures
                        Shear = dees*chi_s ; % Compute shear forces
                    end
                end
                Element_Forces(i,:)=[Moment' Shear'];

            end

            [Row,Column] = size(Element_Forces);
            zeros_row = zeros(1,Column);

            clear Row Column

            assignin('base','Element_Forces',Element_Forces)
            Element_Forces = [zeros_row;Element_Forces];


            W = DISP(:,1);
            assignin('base','W',W)


            for k = 1:analysisObject.STRUCTURE.Number_of_Nodes
                mx = 0. ; my = 0.; mxy = 0.; qx = 0.; qy = 0.;
                ne = 0;
                for iel = 1:analysisObject.STRUCTURE.Number_of_Elements
                    for jel=1:analysisObject.STRUCTURE.number_of_nodes_per_element
                        if analysisObject.connec(iel,jel) == k
                            ne=ne+1;
                            mx = mx + Element_Forces(iel,1);
                            my = my + Element_Forces(iel,2);
                            mxy = mxy + Element_Forces(iel,3);
                            qx = qx + Element_Forces(iel,4);
                            qy = qy + Element_Forces(iel,5);
                        end
                    end
                end
            end
            MX(k,1) = mx/ne;
            MY(k,1) = my/ne;
            MXY(k,1) = mxy/ne;
            QX(k,1) = qx/ne;
            QY(k,1) = qy/ne;

            Results = [MX, MY, MXY, QX, QY];
            figure;
            plot(MX,QX);

            assignin('base','Results',Results)

            %
            figure;
            cmin = min(W);
            cmax = max(W);
            %       caxis([cmin cmax]);
            patch('Faces',analysisObject.connec, 'Vertices', analysisObject.geom, 'FaceVertexCData',W,...
                'Facecolor','interp','Marker','.');
            colorbar;
            toc
        end
    end
end
