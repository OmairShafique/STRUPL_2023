%% This program uses an 4-noded quadrilateral element for the linear
% elastic static analysis of a thick plate in bending
clc %clearing matlab command screen
clear %clearing workspace
tic %starts the clock
% ticBytes(gcp)

%% ----------------------------------------------------------------------
%            ----------  START OF ESSENTIAL INPUT DATA -----------
%------------------------------------------------------------------------

%% Taking essential Inputs
    thickness_of_Plate = input('Please provide Thickness of the plate?: '); %Thickness of the plate

    Elastic_Modulus = input('Please provide Elastic Modulus: '); %Elastic Modulus in

    Poisson_Ratio = input('Please provide Poissons Ratio: ');

%% Declaring Important Global (Strucutral) Variables => i.e. That apply
% to the entire strucutre not to local elements

    % % This Part will ask for File Inputs necessary for further analysis
    % ------------------------------------------------------------------------

        % Loading Coordinates File to
        % Workspace-----------------------------------
            input('Please press enter to find and select your Nodal Coordinates Text File using file explorer:');
            [FileName, Path] = uigetfile('*.txt'); %Will display File Explorer to Choose a text file

            Complete_path = strcat(Path,'\',FileName); %Concatenates string to use in declaring variable to workspace

            %geom
            nodal_coordinate_values = load(Complete_path); %Declaring Variable to workspace

            clear Complete_path FileName Path;%Clearing unnecessary variables


        % Loading Nodal Connectivity for Each
        % Element-----------------------------
            input('Please press enter to find and select your Nodal Connectivity Text File using file explorer:');
            [FileName, Path] = uigetfile('*.txt'); %Will display File Explorer to Choose a text file

            Complete_path = strcat(Path,'\',FileName); %Concatenates string to use in declaring variable to workspace

            %connec
            nodal_connectivity_values = load(Complete_path); %Declaring Variable to workspace

            clear Complete_path FileName Path;%Clearing unnecessary variables
    % ------------------------------------------------------------------------



    % This part will ask for Essential Inputs
%     % ------------------------------------------------------------------------
%         % This input affirms where the environement will be in 2D or 3D
%             %We're not sure what this means exactly,
            dim = input('Would you like to work in 2 or 3 dimensions? Please reply in 2 or 3: ');
% 
%         % This input will ask for the number of nodes per element
%         %Change it to infered values from length and breathe and step
%         %between nodes
% %             number_of_nodes_per_element = input('Please enter the number of nodes per element: ');
% 
%         %This input will ask for the length of the element in the
%         %x-direction
%             Length_of_element = input('Please enter the length of the quadilateral element i.e. in the x-direction: ');
%             
%         %This input will ask for the width of the element in the
%         %y-direction
%             Width_of_element = input('Please enter the width of the quadilateral element i.e. in the y-direction: ');
%             
%             
%         % This input will ask for the number of nodes in the x-direction
%             NXE = input('Please enter the number of nodes in the x-direction: ');
%             
%         % This input will ask for the number of nodes in the x-direction
%             NYE = input('Please enter the number of nodes in the y-direction: '); 
% 
          % This input will ask for the element type
            Element_type = input('Would you like to work with triangle_3, quadratic_4,quadratic_8 or brick_8 ? Please reply in 3,4 or 8: ');
%             
%         % This input will ask for the number of degree of freedoms per node
%             %CHANGE number_of_dof_per_element to number_of_dof_per_node
            number_of_dof_per_node = input('Please enter the number of degrees of freedom per node: ');
            number_of_gauss_point_bending = input('Please enter the number of gauss point bending: ');%ngpb = 3;  number of Gauss points bending
            number_of_gauss_point_shear = input('Please enter the number of gauss point shear: ');%ngps = 2;  number of Gauss points for shear

    % ------------------------------------------------------------------------


    % This part will calculate infered Essential variables
    % ------------------------------------------------------------------------
        
            number_of_nodes_per_element = length(nodal_connectivity_values(1,:)); %Number of Nodes per element
    
            Number_of_Elements = length(nodal_connectivity_values); % Infered Number of Elements

            Number_of_Nodes = length(nodal_coordinate_values); % Infered Number of Nodes in the System
            
            Total_System_Degrees_of_Freedom = Number_of_Nodes * number_of_dof_per_node; %Total System Degrees of Freedom
            
            Degrees_of_Freedom_Per_Element = number_of_nodes_per_element * number_of_dof_per_node; % degrees of freedom per element

            
%             %This variable will calculate the element size in the
%             %x-direction
%             dhx = Length_of_element/NXE;
%             
%             %This variable will calculate the element size in the
%             %x-direction
%             dhy = Length_of_element/NYE;
%             
%             X_origin = 0. ; % x origin of the global coordinate system
%             Y_origin = 0. ; % y origin of the global coordinate system
%% Plotting Mesh Between Nodal Coordinates and Nodal Connectivity
    PlotMesh(nodal_coordinate_values,nodal_connectivity_values);


%% This Part will ask for NON_GLOBAL Input DATA => That applies to local individual elements

% % This Part will ask for File Inputs necessary for further analysis
    % ------------------------------------------------------------------------

        % Loading Boundry Conditions File to Workspace-------------------------------
            input('Please press enter to find and select your Boundry Conditions Text File using file explorer:');
            [FileName, Path] = uigetfile('*.txt'); %Will display File Explorer to Choose a text file

            Complete_path = strcat(Path,'\',FileName); %Concatenates string to use in declaring variable to workspace

            Boundary_Conditions = load(Complete_path); %Declaring Variable to workspace

            clear Complete_path FileName Path;%Clearing unnecessary variables
            
            % Declaring Important Variables from File
                Boundary_Conditions_Degree_of_Freedom = Boundary_Conditions(:,1);%declaring the degree of Freedom of the Boundry Conditions
                    
                    %%TEMPORARY SPACE Populating the NF matrix-------------
                        nf = ones(Number_of_Nodes,number_of_dof_per_node); %nodal freedom matrix set to ones
        
                        node_number_where_restrained = Boundary_Conditions(:,1); %creating a single column matrix with node numbers where dof are restrained
                        dof_x_displacement = Boundary_Conditions(:,2); %dof release matrix from x-disp column in Boundry condition file
                        dof_y_displacement = Boundary_Conditions(:,3); %dof release matrix from y-disp column in Boundry condition file
                        dof_z_displacement = Boundary_Conditions(:,4); %dof release matrix from z-disp column in Boundry condition file
                        dof_x_rotation = Boundary_Conditions(:,5); %dof release matrix from x-rotation column in Boundry condition file
                        dof_y_rotation = Boundary_Conditions(:,6); %dof release matrix from y-rotation column in Boundry condition file
                        dof_z_rotation = Boundary_Conditions(:,7); %dof release matrix from z-rotation column in Boundry condition file
%                         miss = Boundry_Conditions(:,8); %dof release matrix from last column in Boundry condition file.

                        number_of_active_boundary_conditions = length(node_number_where_restrained); %total number of nodes at which restrained present
                        
                        
                        %meant to populate the nf matrix 
                        for i = 1:number_of_active_boundary_conditions
                                if(dim == 2) %checking if 2D 
                                        nf(node_number_where_restrained(i),1) = dof_x_displacement(i); %Changing nf at x_displacmeent
                                        nf(node_number_where_restrained(i),2) = dof_y_displacement(i);
                                        nf(node_number_where_restrained(i),3) = dof_x_rotation(i);
                                        nf(node_number_where_restrained(i),4) = dof_y_rotation(i);
                                else if(dim ==3) %Checking if 3D
                                            nf(node_number_where_restrained(i),1) = dof_x_displacement(i);
                                            nf(node_number_where_restrained(i),2) = dof_y_displacement(i);
                                            nf(node_number_where_restrained(i),3) = dof_z_displacement(i);
                                            nf(node_number_where_restrained(i),4) = dof_x_rotation(i);
                                            nf(node_number_where_restrained(i),5) = dof_y_rotation(i);
                                            nf(node_number_where_restrained(i),6) = dof_z_rotation(i);
                                            nf(node_number_where_restrained(i),7) = miss(i);
                                    end
                                end                             
                        end
                        
                        
                        %meant to count total number of active degrees of
                        %freedom 
                        total_numbers_of_active_dof = 0;
                        
                        for i = 1:length(nf)
                            for j = 1:3
                                if(nf(i,j) == 1) %checking if nf has 1 or 0
                                    total_numbers_of_active_dof = total_numbers_of_active_dof + 1; %if 1 then increase value by 1
                                    nf(i,j)= total_numbers_of_active_dof;
                                end
                            end
                        end
                         % Counting of the free degrees of freedom



disp ('Nodal freedom')
nf
disp ('Total number of active degrees of freedom')
 total_numbers_of_active_dof
                        
                        clear node_number_where_active
                        clear dof_x_displacement
                        clear dof_y_displacement
                        clear dof_z_displacement
                        clear dof_x_rotation
                        clear dof_y_rotation
                        clear dof_z_rotation
                        clear miss        
        %Loading External Loads file to workspace--------------------------
            input('Please press enter to find and select your External Loads Text File using file explorer:');
            [FileName, Path] = uigetfile('*.txt'); %Will display File Explorer to Choose a text file

            Complete_path = strcat(Path,'\',FileName); %Concatenates string to use in declaring variable to workspace

            External_load = load(Complete_path); %Declaring Variable to workspace

            clear Complete_path FileName Path;%Clearing unnecessary variables
            
            %Declaring Important variables from File
                Node_number_of_external_load= External_load(:,1);
                Load = zeros(Number_of_Nodes,3); 
                Force_x_direction=External_load(:,2);
                Force_y_direction=External_load(:,3);
                Force_z_direction=External_load(:,4);
                Nodal_load=length(Node_number_of_external_load);
%                 
                for i=1:Nodal_load
%               
                        Load(Node_number_of_external_load(i),1)=Force_x_direction(i);
                        Load(Node_number_of_external_load(i),2)=Force_y_direction(i);
                        Load(Node_number_of_external_load(i),3)=Force_z_direction(i);                       
                end 
%                 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% End of input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Assemble the global force vector
%
Global_force_vector=zeros(total_numbers_of_active_dof,1);
for i=1: Number_of_Nodes
    for j=1:number_of_dof_per_node
        if nf(i,j) ~= 0
            Global_force_vector(nf(i,j))= Load(i,j);
        end
    end
end
% 
%
%% %%%%%%%%% Numerical integration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Form the matrix containing the abscissas and the weights of Gauss points
%
sampb=gauss(number_of_gauss_point_bending); samps=gauss(number_of_gauss_point_shear);                                                                                                  
%
% Initialize the global stiffness matrix to zero
Global_stiffness_matrix = zeros(total_numbers_of_active_dof, total_numbers_of_active_dof);

for i=1: Number_of_Elements
coord=zeros(number_of_nodes_per_element,dim);
for k=1: number_of_nodes_per_element
    for j=1:dim
        coord(k,j)=nodal_coordinate_values(nodal_connectivity_values(i,k),j);
    end
end
    %
l=0;
for k=1: number_of_nodes_per_element
   for j=1:number_of_dof_per_node
        l=l+1;
        g(l)=nf(nodal_connectivity_values(i,k),j);
    end
end
end
% [coord,g] = platelem_q8(i) ; % coordinates of the nodes of element i,
% and its steering vector
keb=zeros(Degrees_of_Freedom_Per_Element,Degrees_of_Freedom_Per_Element); % Initialize the element bending
% stiffness matrix to zero
kes=zeros(Degrees_of_Freedom_Per_Element,Degrees_of_Freedom_Per_Element); % Initialize the element Shear
% stiffness matrix to zero
% Integrate element bending stiffness and assemble it in global matrix
%
for ig=1: number_of_gauss_point_bending
wi = sampb(ig,2);
for jg=1: number_of_gauss_point_bending
wj=sampb(jg,2);
[der,fun] = fmquad(sampb,Element_type,dim, ig,jg); % Derivative of shape functions
% in local coordinates
jac=der.*coord; % Compute Jacobian matrix
d=det(jac); % Compute the determinant of
% Jacobian matrix
jac1=inv(jac); % Compute inverse of the Jacobian
deriv=jac1*der; % Derivative of shape functions
% in global coordinates
beeb=formbeeb(deriv,number_of_nodes_per_element,Degrees_of_Freedom_Per_Element);% Form matrix [B]
deeb=formdeeb(Elastic_Modulus,Poisson_Ratio,thickness_of_Plate);
keb=keb + d*wi*wj*beeb'*deeb*beeb; % Integrate stiffness matrix
end
end
%  Global_stiffness_matrix=form_Global_stiffness_matrix(Global_stiffness_matrix,keb, g); % assemble global stiffness matrix
for iel=1:Degrees_of_Freedom_Per_Element
    if g(iel) ~= 0
        for j=1: Degrees_of_Freedom_Per_Element
            if g(j) ~= 0
                Global_stiffness_matrix(g(iel),g(j))= Global_stiffness_matrix(g(iel),g(j)) + keb(iel,j);
            end
        end
    end

% Integrate element Shear stiffness and assemble it in global matrix
%
for ig=1: number_of_gauss_point_shear
wi = samps(ig,2);
for jg=1: number_of_gauss_point_shear
wj=samps(jg,2);
[der,fun] = fmquad(samps, ig,jg); % Derivative of shape functions
% in local coordinates
jac=der*coord; % Compute Jacobian matrix
d=det(jac); % Compute determinant of
            % Jacobian matrix
jac1=inv(jac); % Compute inverse of the  Jacobian

deriv=jac1*der; % Derivative of shape functionsin global coordinates

dees=formdees(Elastic_Modulus,Poisson_Ratio,thickness_of_Plate);
bees=formbees(deriv,fun,number_of_nodes_per_element,Degrees_of_Freedom_Per_Element); % Form matrix [B]
kes=kes+(5/6)*d*wi*wj*bees'*dees*bees; % Integrate stiffness matrix
end
end
for iel=1:Degrees_of_Freedom_Per_Element
    if g(iel) ~= 0
        for j=1: Degrees_of_Freedom_Per_Element
            if g(j) ~= 0
                Global_stiffness_matrix(g(iel),g(j))= Global_stiffness_matrix(g(iel),g(j)) + kes(iel,j); 
            end
        end
    end
end
end
%
%
%% %%%%%%%%%%%%%%%%%%%%% End of assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
delta = Global_stiffness_matrix\Global_force_vector % solve for unknown displacements
format short e
disp('node w_disp x_slope y_slope ') %
%
for i=1: Number_of_Nodes %
if nf(i,1) == 0 %
w_disp =0.; %
else
w_disp = delta(nf(i,1)); %
end
%
if nf(i,2) == 0 %
x_slope = 0.; %
else
x_slope = delta(nf(i,2)); %
end
%
if nf(i,3) == 0 %
y_slope = 0.; %
else
y_slope = delta(nf(i,3)); %
end
disp([i w_disp x_slope y_slope]) % Display displacements of each node
DISP(i,:) = [ w_disp x_slope y_slope];
end
%
%
%
ngp=1; % Calculate moments and shear forces
% the center of each element
samp=gauss(ngp);
%
% for i=1:Number_of_Elements
% [coord,g] = platelem_q8(i); % coordinates of the nodes of element i,
% % and its steering vector
for i=1: Number_of_Elements
coord=zeros(number_of_nodes_per_element,dim);
for k=1: number_of_nodes_per_element
    for j=1:2
        coord(k,j)=nodal_coordinate_values(nodal_connectivity_values(i,k),j);
    end
end
    %
l=0;
for k=1: number_of_nodes_per_element
    for j=1:number_of_dof_per_node
        l=l+1;
        g(l)=nf(nodal_connectivity_values(i,k),j);
    end
end
element_displacement=zeros(Degrees_of_Freedom_Per_Element,1); % Initialize element displacement to zero
for m=1:Degrees_of_Freedom_Per_Element %
if g(m)==0 %
element_displacement(m)=0.; %
else %
element_displacement(m)=delta(g(m)); % Retrieve element displacement from the
                                     % global displacement vector
end
end
%
for ig=1: ngp
wi = samp(ig,2);
for jg=1: ngp
wj=samp(jg,2);
[der,fun] = fmquad(samp, ig,jg); % Derivative of shape functions
                                 % in local coordinates
jac=der*coord;                   % Compute Jacobian matrix
d=det(jac);                      % Compute the determinant of
                                 % Jacobian matrix
jac1=inv(jac);                   % Compute inverse of the Jacobian
deriv=jac1*der;                  % Derivative of shape functions
                                 % in global coordinates
%
beeb=formbeeb(deriv,number_of_nodes_per_element,Degrees_of_Freedom_Per_Element);% Form matrix [B]
chi_b = beeb*element_displacement ; % compute bending curvatures
Moment = deeb*chi_b ;               % Compute moments
bees=formbees(deriv,fun,number_of_nodes_per_element,Degrees_of_Freedom_Per_Element); % Form matrix [B]
chi_s = bees*element_displacement ; % compute shear curvatures
Shear = dees*chi_s ;                % Compute shear forces
end
end
Element_Forces(i,:)=[Moment' Shear'];
end
%
W = DISP(:,1);
[MX, MY, MXY, QX, QY] = Forces_at_nodes_plate(Element_Forces);
%

cmin = min(W);
cmax = max(W);
caxis([cmin cmax]);
patch('Faces',nodal_connectivity_values, 'Vertices', nodal_coordinate_values, 'FaceVertexCData',W,...
'Facecolor','interp','Marker','.');
colorbar;
toc
tocBytes(gcp)






