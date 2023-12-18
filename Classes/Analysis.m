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
        ngps = 0; % Number of Gauss Points Sheer
        ngpb = 0; % Number of Gauss Points Bending

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
        function analysis = Analysis(thickness_of_plate,number_of_dof_per_node,dim,Boundary_Conditions,External_Load,nodal_coordinate_values,nodal_connectivity_values)
            analysis.thickness_of_plate = thickness_of_plate;
            analysis.number_of_dof_per_node = number_of_dof_per_node;
            analysis.dim = dim;
            analysis.Boundary_Conditions = Boundary_Conditions;
            analysis.External_Load = External_Load;

            analysis.STRUCTURE = Structure(nodal_coordinate_values,nodal_connectivity_values);
        end

        % Main Class where all logic will follow
        function Engine()

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

        % Populating Load Matric
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

        % GaussFunction
        function [samp]=gauss(ngp)
            %
            % This function returns the abscissas and weights of the Gauss
            % points for ngp equal up to 4
            %
            %
            samp=zeros(ngp,2);
            %
            if ngp==1
                samp=[0. 2];
            elseif ngp==2
                samp=[-1./sqrt(3) 1.;...
                    1./sqrt(3) 1.];
            elseif ngp==3
                samp= [-.2*sqrt(15.) 5./9; ...
                    0. 8./9.;...
                    .2*sqrt(15.) 5./9];
            elseif ngp==4
                samp= [-0.861136311594053 0.347854845137454; ...
                    -0.339981043584856 0.652145154862546; ...
                    0.339981043584856 0.652145154862546; ...
                    0.861136311594053 0.347854845137454];
            end
            %
            % End function Gauss
        end

        % Fg_gravity_Matrix_Calculator
        function fg_gravity_Final_end = fg_matrix_calculator()

            %% Finding a matrix that will provide the number of times a node is repeated and work with fg_matrix to
            % find the gravity loading matrix in tune witht the Global_Force_Vector so
            % these two can be added.

            [rows,columns] = size(this.nodal_connectivity_values);
            index_length = rows*columns;


            if this.STRUCTURE.Element_Type == 3
                this.STRUCTURE.number_of_nodes_per_element = 3;
            else
                this.STRUCTURE.number_of_nodes_per_element = length(this.nodal_connectivity_values(1,:)); %Number of Nodes per element
            end

            fg_gravity = zeros(this.STRUCTURE.number_of_elements * this.STRUCTURE.number_of_nodes_per_element * this.number_of_dof_per_node,1);


            Gravity_Load = [0 ,-this.gamma]';

            current_row = 1;
            for iel = 1:number_of_elements

                [bee,fun_3,g,A,d_3] = elem_T3(iel); % HAVE TO MAKE A FUNCTION FOR THIS TOO

                fg_gravity(current_row: current_row + this.STRUCTURE.number_of_dof_per_node*number_of_nodes_per_element - 1) = ...
                    (fun_3 * Gravity_Load * d_3) * this.thickness_of_plate * (-1/3);

                current_row = current_row + this.STRUCTURE.number_of_dof_per_node * this.STRUCTURE.number_of_nodes_per_element;
            end

            [rows,columns] = size(fg_gravity);

            fg_gravity_x2 = zeros(rows/2,3);

            index = 1;

            for i = 1:rows/2
                for j = 2:3
                    fg_gravity_x2(i,j) = fg_gravity(index);
                    index = index + 1;
                end
            end

            [rows,columns] = size(this.nodal_connectivity_values);
            index = 1;

            for i = 1:rows
                for j = 1:columns
                    fg_gravity_x2(index,1) = this.nodal_connectivity_values(i,j);
                    index = index + 1;
                end
            end

            fg_gravity_Final = zeros(Number_of_Nodes,3);
            [rows,columns] = size(fg_gravity_x2);
            index = 1;

            for i = 1:this.STRUCTURE.Number_of_Nodes
                for j = 1:rows
                    if fg_gravity_x2(j,1) == i
                        fg_gravity_Final(i,1) = i;
                        fg_gravity_Final(i,2:3) = fg_gravity_Final(i,2:3) + fg_gravity_x2(j,2:3);
                    end
                end
            end

            [rows,columns] = size(fg_gravity_Final);
            index = 1;
            fg_gravity_Final_end = zeros(rows*(columns-1),1);

            for i = 1:rows
                for j = 2:columns
                    fg_gravity_Final_end(index) = fg_gravity_Final(i,j);
                    index = index + 1;
                end
            end


        end

        % form deeb
        function [deeb] = formdeeb(Elatic_Modulus,Poisson_Ratio,thickness_of_Plate)
            DR= Elatic_Modulus*(thickness_of_Plate^3)/(12*(1.-Poisson_Ratio*Poisson_Ratio));
            deeb=DR*[1 Poisson_Ratio 0. ;...
                Poisson_Ratio 1 0. ;...
                0. 0. (1.-Poisson_Ratio)/2] ;
        end

        % form dees
        function [dees] = formdees(Elatic_Modulus,Poisson_Ratio,thickness_of_Plate)
            G= Elatic_Modulus/(2*(1.+Poisson_Ratio));
            dees=G* [thickness_of_Plate     0 ;...
                0     thickness_of_Plate];
        end

        % elem_T3
        function [bee,fun_3,g,A,d_3] = elem_T3(iel)
            %
            % This function returns the coordinates of the nodes of element i
            % and its steering vector
            %

            % if Element_Type==3
            x1 = this.geom(this.connec(iel,1),1); y1 = this.geom(this.connec(iel,1),2);
            x2 = this.geom(this.connec(iel,2),1); y2 = this.geom(this.connec(iel,2),2);
            x3 = this.geom(this.connec(iel,3),1); y3 = this.geom(this.connec(iel,3),2);
            %
            %
            d_3=det([1 x1 y1; ...
                1 x2 y2; ...
                1 x3 y3]);

            A = (0.5)*det([1 x1 y1; ...
                1 x2 y2; ...
                1 x3 y3]);


            m11 = (x2*y3 - x3*y2)/(2*A);
            m21 = (x3*y1 - x1*y3)/(2*A);
            m31 = (x1*y2 - y1*x2)/(2*A);
            m12 = (y2 - y3)/(2*A);
            m22 = (y3 - y1)/(2*A);
            m32 = (y1 - y2)/(2*A);
            m13 = (x3 - x2)/(2*A);
            m23 = (x1 - x3)/(2*A);
            m33 = (x2 -x1)/(2*A);

            bee = [ m12 0 m22 0 m32 0; ...
                0 m13 0 m23 0 m33; ...
                m13 m12 m23 m22 m33 m32] ;

            fun_3= [(m11+m12+m13)  0;...
                0  (m11+m12+m13);...
                (m21+m22+m23)  0;...
                0  (m21+m22+m23);...
                (m31+m32+m33)  0;...
                0  (m31+m32+m33)];
            %

            l=0;
            for k=1: number_of_nodes_per_element
                for j=1:this.STRUCTURE.number_of_dof_per_node
                    l=l+1;
                    g(l)=nf_g(this.connec(iel,k),j);
                end
            end



        end

        % Plate_q4
        function [coord,g] = platelem_q4(i)
            %
            % This function returns the coordinates of the nodes of element i
            % and its steering vector
            % %

            % To cchange the size of the problem or change the elastic properties
            % ALTER the PlateQ8_input_module.m
            %
            %in form_KK where the KK() is trying to access a g() with a dimension = eldof but
            %is restricted to only number_of_nodes_per_element or dim

            coord=zeros(this.STRUCTURE.number_of_nodes_per_element,this.dim);
            for k=1: number_of_nodes_per_element
                for j=1:this.dim
                    coord(k,j)=this.geom(this.connec(i,k),j);
                end
            end
            %
            l=0;
            g = 0;
            for k=1:number_of_nodes_per_element
                for j=1:this.STRUCTURE.number_of_dof_per_node
                    l=l+1;
                    g(l)=nf_g(this.connec(i,k),j);
                end
            end

        end

        % fmquad
        function [der,fun] = fmquad(samp, ig,jg,iel)%,lg

            %% To be called only for Element_Type == 8 or == 4
            %


            % This function
            % returns the vector of the shape function and their
            % derivatives with respect to xi and eta at the gauss points for
            % an 8-nodded quadrilateral
            %
            % Element_Type=8;
            % Element_Type=3;
            % Element_Type=4;
            % dim=2;
            % dim=3;
            xi=samp(ig,1);
            eta=samp(jg,1);
            % zeta=samp(lg,1);
            etam=(1.-eta);
            etap=(1.+eta);
            xim=(1.-xi);
            xip=(1.+xi);
            % zetam=(1.-zeta);
            % zetap=(1.+zeta);

            %

            switch (Element_Type)

                case (3) % for triangular element
                    if this.Element_Type==3 && this.dim==2 && this.ngpb ~=0 && this.ngps ~=0
                        % shape functions
                        fun(1) = 1. - xi - eta;
                        fun(2) =  xi;
                        fun(3) =  eta;
                        % derivatives


                        der(1,1)= -1.; der(1,2)=-1;
                        der(2,1)= 1.; der(2,2)= 0;
                        der(3,1)=0; der(3,2)=1;
                    else
                        x1 = this.geom(this.connec(iel,1),1); y1 = this.geom(this.connec(iel,1),2);
                        x2 = this.geom(this.connec(iel,2),1); y2 = this.geom(this.connec(iel,2),2);
                        x3 = this.geom(this.connec(iel,3),1); y3 = this.geom(this.connec(iel,3),2);

                        A = (0.5)*det([1 x1 y1; ...
                            1 x2 y2; ...
                            1 x3 y3]);

                        m11 = (x2*y3 - x3*y2)/(2*A); % reanme to ders
                        m21 = (x3*y1 - x1*y3)/(2*A);
                        m31 = (x1*y2 - y1*x2)/(2*A);
                        m12 = (y2 - y3)/(2*A);
                        m22 = (y3 - y1)/(2*A);
                        m32 = (y1 - y2)/(2*A);
                        m13 = (x3 - x2)/(2*A);
                        m23 = (x1 - x3)/(2*A);
                        m33 = (x2 -x1)/(2*A);

                        % der(1,1)  = m11;
                        % der(2,1)  = m21;
                        % der(3,1)  = m31;
                        der(1,1)  = m12;
                        der(1,2)  = m22;
                        der(2,1)  = m32;
                        der(2,2)  = m13;
                        der(3,1)  = m23;
                        der(3,2)  = m33;

                        fun = [(m11+m12+m13)  0;...
                            0  (m11+m12+m13);...
                            (m21+m22+m23)  0;...
                            0  (m21+m22+m23);...
                            (m31+m32+m33)  0;...
                            0  (m31+m32+m33)];

                    end
                case (4) % for rectangular element

                    % shape functions
                    fun(1)=0.25*xim*etam;
                    fun(2)=0.25*xip*etam;
                    fun(3)=0.25*xip*etap;
                    fun(4)=0.25*xim*etap;
                    % derivatives
                    der(1,1)=-0.25*etam; der(1,2)=-0.25*xim;
                    der(2,1)=0.25*etam;  der(2,2)=-0.25*xip;
                    der(3,1)=0.25*etap;  der(3,2)=0.25*xip;
                    der(4,1)=-0.25*etap; der(4,2)=0.25*xim;

                case (8) % for 8-noded Element

                    if this.Element_Type==8 && this.dim==2
                        % shape functions
                        fun(1) = -0.25*xim*etam*(1.+ xi + eta);
                        fun(2) = 0.5*(1.- xi^2)*etam;
                        fun(3) = -0.25*xip*etam*(1. - xi + eta);
                        fun(4) = 0.5*xip*(1. - eta^2);
                        fun(5) = -0.25*xip*etap*(1. - xi - eta);
                        fun(6) = 0.5*(1. - xi^2)*etap;
                        fun(7) = -0.25*xim*etap*(1. + xi - eta);
                        fun(8) = 0.5*xim*(1. - eta^2);
                        % derivatives
                        der(1,1)=0.25*etam*(2.*xi + eta); der(1,2)=-1.*etam*xi;
                        der(1,3)=0.25*etam*(2.*xi-eta); der(1,4)=0.5*(1-eta^2);
                        der(1,5)=0.25*etap*(2.*xi+eta); der(1,6)=-1.*etap*xi;
                        der(1,7)=0.25*etap*(2.*xi-eta); der(1,8)=-0.5*(1.-eta^2);
                        %
                        der(2,1)=0.25*xim*(2.*eta+xi); der(2,2)=-0.5*(1. - xi^2);
                        der(2,3)=-0.25*xip*(xi-2.*eta); der(2,4)=-1.*xip*eta;
                        der(2,5)=0.25*xip*(xi+2.*eta); der(2,6)=0.5*(1.-xi^2);
                        der(2,7)=-0.25*xim*(xi-2.*eta); der(2,8)=-1.*xim*eta;
                    else
                        if this.Element_Type==8 && this.dim==3  %brick element
                            % shape functions
                            fun(1) = 0.125*xim*etam*zetam;
                            fun(2) = 0.125*xip*etam*zetam;
                            fun(3) = 0.125*xip*etap*zetam;
                            fun(4) = 0.125*xim*etap*zetam;
                            fun(5) = 0.125*xim*etam*zetap;
                            fun(6) = 0.125*xip*etam*zetap;
                            fun(7) = 0.125*xip*etap*zetap;
                            fun(8) = 0.125*xim*etap*zetap;
                            % derivatives
                            der(1,1)=0.125*etam*zetam; der(1,2)=0.125*etam*zetam;
                            der(1,3)=0.125*etap*zetam; der(1,4)=0.125*etap*zetam;
                            der(1,5)=0.125*etam*zetap; der(1,6)=0.125*etam*zetap;
                            der(1,7)=0.125*etap*zetap; der(1,8)=0.125*etap*zetap;
                            %
                            der(2,1)=0.125*xim*zetam; der(2,2)=0.125*xip*zetam;
                            der(2,3)=0.125*xip*zetam; der(2,4)=0.125*xim*zetam;
                            der(2,5)=0.125*xim*zetap; der(2,6)=0.125*xip*zetap;
                            der(2,7)=0.125*xip*zetap; der(2,8)=0.125*xim*zetap;
                            %
                            der(3,1)=0.125*xim; der(3,2)=0.125*xip;
                            der(3,3)=0.125*xip; der(3,4)=0.125*xim;
                            der(3,5)=0.125*xim; der(3,6)=0.125*xip;
                            der(3,7)=0.125*xip; der(3,8)=0.125*xim;
                        end
                    end
            end
        end

        % form beeb and
        function [beeb] = formbeeb(deriv)
            %
            % This function assembles the matrix [beeb] from the
            % derivatives of the shape functions in global coordinates
            % for a thick plate element (bending action)
            %
            beeb=zeros(3,this.Degrees_of_Freedom_Per_Element);
            for m=1:this.number_of_nodes_per_element
                k=this.number_of_dof_per_node*m;
                j=k-1;
                x=deriv(1,m);
                beeb(1,j)=x;
                beeb(3,k)=x;
                y=deriv(2,m);
                beeb(2,k)=y;
                beeb(3,j)=y;
            end
        end

        % form bees and
        function [bees] = formbees(deriv,fun,Degrees_of_Freedom_Per_Element)

            bees=zeros(2,Degrees_of_Freedom_Per_Element);
            for m=1:this.number_of_nodes_per_element
                k=this.number_of_dof_per_node*m;
                j=k-1;
                i=k-2;
                x=deriv(1,m); y=deriv(2,m);
                bees(2,i)=-x;
                bees(1,i)=-y;
                bees(1,k) = fun(m);
                bees(2,j) = fun(m);
            end


        end

        % form KK Matrix
        function [Global_stiffness_matrix] = form_KK(Global_stiffness_matrix, kg, g)
            % This function assembles the global stiffness matrix
            % Degrees_of_Freedom_Per_Elementdof = 12;
            % This function assembles the global stiffness matrix

            for i=1:Degrees_of_Freedom_Per_Element
                if g(i) ~= 0
                    for j=1: Degrees_of_Freedom_Per_Element
                        if g(j) ~= 0
                            Global_stiffness_matrix(g(i),g(j))= Global_stiffness_matrix(g(i),g(j)) + kg(i,j);
                        end
                    end
                end
            end
        end

        % Plot Mesh
        function PlotMesh(coordinates,nodes)

            %%INTRODUCE 3D
            %--------------------------------------------------------------------------
            % Code written by : Siva Srinivas Kolukula                                |
            %                   Senior Research Fellow                                |
            %                   Structural Mechanics Laboratory                       |
            %                   Indira Gandhi Center for Atomic Research              |
            %                   India                                                 |
            % E-mail : allwayzitzme@gmail.com                                         |
            %          http://sites.google.com/site/kolukulasivasrinivas/             |
            %--------------------------------------------------------------------------
            %--------------------------------------------------------------------------
            % Purpose:
            %         To plot the Finite Element Method Mesh
            % Synopsis :
            %           PlotMesh(coordinates,nodes)
            % Variable Description:
            %           coordinates - The nodal coordinates of the mesh
            %           -----> coordinates = [node X Y]
            %           nodes - The nodal connectivity of the elements
            %           -----> nodes = [node1 node2......]
            %--------------------------------------------------------------------------

            nel = length(nodes) ;                  % number of elements
            nnd = length(coordinates) ;          % total number of nodes in system
            nne = size(nodes,2);                % number of nodes per element
            %
            % Initialization of the required matrices
            X = zeros(nne,nel) ;
            Y = zeros(nne,nel) ;

            for iel=1:nel
                for i=1:nne
                    nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
                    X(i,iel)=coordinates(nd(i),1);    % extract x value of the node
                    Y(i,iel)=coordinates(nd(i),2);    % extract y value of the node
                end
            end
            if nne==8

                % Plotting the FEM mesh, diaplay Node numbers and Element numbers
                f1 = figure ;
                set(f1,'name','Mesh','numbertitle','off') ;
                plot(X,Y,'k')
                fill(X,Y,'w')

                title('Finite Element Mesh');
                axis off ;
                k = nodes(:,1:end);
                nd = k' ;
                for i = 1:nel
                    text(X(:,i),Y(:,i),int2str(nd(:,i)),'fontsize',8,'color','k');
                    text(sum(X(:,i))/8,sum(Y(:,i))/8,int2str(i),'fontsize',10,'color','r') ;
                end
            elseif nne==4
                % Plotting the FEM mesh, diaplay Node numbers and Element numbers
                f1 = figure ;
                set(f1,'name','Mesh','numbertitle','off') ;
                plot(X,Y,'k')
                fill(X,Y,'w')

                title('Finite Element Mesh') ;
                axis off ;
                k = nodes(:,1:end);
                nd = k' ;
                for i = 1:nel
                    text(X(:,i),Y(:,i),int2str(nd(:,i)),'fontsize',8,'color','k');
                    text(sum(X(:,i))/4,sum(Y(:,i))/4,int2str(i),'fontsize',10,'color','r') ;
                end
            elseif nne==3
                % Plotting the FEM mesh, diaplay Node numbers and Element numbers
                f1 = figure ;
                set(f1,'name','Mesh','numbertitle','off') ;
                plot(X,Y,'k')
                fill(X,Y,'w')

                title('Finite Element Mesh') ;
                axis off ;
                k = nodes(:,1:end);
                nd = k' ;
                for i = 1:nel
                    text(X(:,i),Y(:,i),int2str(nd(:,i)),'fontsize',8,'color','k');
                    text(sum(X(:,i))/3,sum(Y(:,i))/3,int2str(i),'fontsize',10,'color','r') ;
                end
            end

        end

        % Elastic Solve
        function [delta,Reaction,node_displacement,SIGMA,STRAIN] = Elastic_solve()


            % dealta is related to nodal dislplacement
            delta = KK\fg ; % solve for unknown displacements
            Reaction=KK*delta-fg;
            %
            node_displacement=zeros(nnd,2);
            %
            for i=1: nnd %
                if nf(i,1) == 0 %
                    x_disp =0.; %
                else
                    x_disp = delta(nf(i,1)); %
                end
                %
                if nf(i,2) == 0 %
                    y_disp = 0.; %
                else
                    y_disp = delta(nf(i,2)); %
                end
                node_displacement(i,:) =[x_disp y_disp];
            end
            %
            % Retrieve the x_coord and y_disp of the nodes located on the neutral axis
            %
            k = 0;
            for i=1:nnd
                if geom(i,2)== 0.
                    k=k+1;
                    x_coord(k) = geom(i,1);
                    vertical_disp(k)=node_displacement(i,2);
                end
            end
            %% Compute stress and strain %%
            %
            for i=1:nel
                [bee,g,A] = elem_T3(i); % Form strain matrix, and steering vector
                eld=zeros(eldof,1); % Initialize element displacement to zero
                for m=1:eldof
                    if g(m)==0
                        eld(m)=0.;
                    else %
                        eld(m)=delta(g(m)); % Retrieve element displacement
                    end
                end
                %
                eps=bee*eld; % Compute strains
                STRAIN(i,:)=eps ; % Store strains for all elements
                sigma=dee*eps; % Compute stresses
                SIGMA(i,:)=sigma ; % Store stress for all elements
            end

        end

        % Form B Matrix
        function [B_Matrix] = form_B_Matrix()
            % This function assembles the matrix B for plane problems from the
            % derivatives of the shape functions in global coordinates
            % for a thick plate element (bending action)

            if this.STRUCTURE.Element_type==3 && this.ngpb==0

                x1 = this.geom(this.connec(iel,1),1); y1 = this.geom(this.connec(iel,1),2);
                x2 = this.geom(this.connec(iel,2),1); y2 = this.geom(this.connec(iel,2),2);
                x3 = this.geom(this.connec(iel,3),1); y3 = this.geom(this.connec(iel,3),2);
                %
                A = (0.5)*det([1 x1 y1; ...
                    1 x2 y2; ...
                    1 x3 y3]);
                %
                m11 = (x2*y3 - x3*y2)/(2*A);
                m21 = (x3*y1 - x1*y3)/(2*A);
                m31 = (x1*y2 - y1*x2)/(2*A);
                m12 = (y2 - y3)/(2*A);
                m22 = (y3 - y1)/(2*A);
                m32 = (y1 - y2)/(2*A);
                m13 = (x3 - x2)/(2*A);
                m23 = (x1 - x3)/(2*A);
                m33 = (x2 -x1)/(2*A);
                %
                B_Matrix = [ m12 0 m22 0 m32 0; ...
                    0 m13 0 m23 0 m33; ...
                    m13 m12 m23 m22 m33 m32] ;
            end
            if this.Element_type~=3 && this.ngpb~=0

                B_Matrix=zeros(3,this.STRUCTURE.Degrees_of_Freedom_Per_Element);
                for m=1:this.number_of_nodes_per_element
                    k=3*m;
                    j=k-1;
                    x=deriv(1,m);
                    B_Matrix(1,j)=x;
                    B_Matrix(3,k)=x;
                    y=deriv(2,m);
                    B_Matrix(2,k)=y;
                    B_Matrix(3,j)=y;
                end
            end
            if this.Element_type~=3 && this.ngps~=0
                B_Matrix=zeros(2,this.STRUCTURE.Degrees_of_Freedom_Per_Element);
                for m=1:this.number_of_nodes_per_element
                    k=3*m;
                    j=k-1;
                    i=k-2;
                    x=deriv(1,m); y=deriv(2,m);
                    B_Matrix(2,i)=-x;
                    B_Matrix(1,i)=-y;
                    B_Matrix(1,k) = fun(m);
                    B_Matrix(2,j) = fun(m);
                end
            end
            % End function form_B_Matrix



        end

        % Form D Matrix
        function [D_Matrix] = form_D_Matrix ()
            %
            % This function forms the elasticity matrix in a plate element
            % concerning also shear and bending actions
            %
            if this.Element_type==3 && this.ngpb==0
                c = this.Elastic_Modulus/(1.- this.Poisson_Ratio * this.Poisson_Ratio);
                %
                D_Matrix=c*[1 this.Poisson_Ratio 0. ;...
                    this.Poisson_Ratio 1 0. ;...
                    0. 0. .5*(1.- this.Poisson_Ratio)];

            elseif this.Element_type~=3 && this.ngpb~=0
                DR= this.Elastic_Modulus*(this.thickness_of_Plate^3)/(12*(1.-this.Poisson_Ratio*this.Poisson_Ratio));
                %
                D_Matrix=DR*[1 this.Poisson_Ratio 0. ;...
                    this.Poisson_Ratio 1 0. ;...
                    0. 0. (1.- this.Poisson_Ratio)/2] ;
            elseif this.Element_type~=3 && this.ngps~=0
                G= this.Elatic_Modulus/(2*(1. + this.Poisson_Ratio));
                %
                D_Matrix=G* [this.thickness_of_Plate     0 ;...
                    0     this.thickness_of_Plate];
            end

        end

        % Forces at nodes plates
        function [MX, MY, MXY, QX, QY]=Forces_at_nodes_plate(Element_Forces)
            % This function averages the stresses at the nodes
            for k = 1:this.STRUCTURE.Number_of_Nodes
                mx = 0. ; my = 0.; mxy = 0.; qx = 0.; qy = 0.;
                ne = 0;
                for iel = 1:this.STRUCTURE.Number_of_Elements
                    for jel=1:this.STRUCTURE.number_of_nodes_per_element
                        if this.connec(iel,jel) == k
                            ne=ne+1;
                            mx = mx + Element_Forces(iel,1);
                            my = my + Element_Forces(iel,2);
                            mxy = mxy + Element_Forces(iel,3);
                            qx = qx + Element_Forces(iel,4);
                            qy = qy + Element_Forces(iel,5);
                        end
                    end
                end
                MX(k,1) = mx/ne;
                MY(k,1) = my/ne;
                MXY(k,1) = mxy/ne;
                QX(k,1) = qx/ne;
                QY(k,1) = qy/ne;
            end
        end
    
        % 
    
    end
end