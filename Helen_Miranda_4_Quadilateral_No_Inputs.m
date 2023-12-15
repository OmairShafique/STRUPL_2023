%% This program uses an 4-noded quadrilateral element for the linear
% elastic statis analysis of a thick plate in bending
clc %clearing matlab command screen
tic %starts the clock

%% ----------------------------------------------------------------------
%            ----------  START OF ESSENTIAL INPUT DATA -----------
%------------------------------------------------------------------------

%% Declaring Important Global (Strucutral) Variables => i.e. That apply
% to the entire strucutre not to local elements

%this line is from Omair and not helen
%this change should show up in omair and not main

%this line is from helens
% This part will calculate infered Essential variables
% ------------------------------------------------------------------------
%comit
global geom ngpb nf_g  dees ngps deeb connec total_numbers_of_active_dof Number_of_Nodes Number_of_Elements Degrees_of_Freedom_Per_Element number_of_nodes_per_element nf dim number_of_dof_per_node

number_of_dof_per_node = evalin('base','nodof');
dim = evalin('base','dim');
Boundary_Conditions = evalin('base','Boundary_Conditions');
Elastic_Modulus = evalin('base','Elastic_Modulus');
External_Load = evalin('base','External_Load');
Length = evalin('base','Length_of_Element');
Width = evalin('base','Width_of_element');
NXE = evalin('base','NXE');
NYE = evalin('base','NYE');
Poissons_Ratio = evalin('base','Possions_Ratio');
thickness_of_plate = evalin('base','thickness_of_plate');
gamma = evalin('base','weigth_density'); %Insert weight gamma
node = evalin('base','nodof');
nodal_coordinate_values = evalin('base','nodal_coordinate_values');
nodal_connectivity_values = evalin('base','nodal_connectivity_values');

geom = evalin('base','nodal_coordinate_values'); %Geom Matrix
connec = evalin('base','nodal_connectivity_values'); %connec Matrix

% This input will ask for the element type

Element_Type = evalin('base','Element_Type');

if Element_Type == 3
    number_of_nodes_per_element = 3;
else
    number_of_nodes_per_element = length(nodal_connectivity_values(1,:)); %Number of Nodes per element
end

assignin('base','number_of_nodes_per_element',number_of_nodes_per_element);

Number_of_Elements = length(nodal_connectivity_values); % Infered Number of Elements
assignin('base','Number_of_Elements',Number_of_Elements);


Number_of_Nodes = length(nodal_coordinate_values); % Infered Number of Nodes in the System
assignin('base','Number_of_Nodes',Number_of_Nodes);



Total_System_Degrees_of_Freedom = Number_of_Nodes * number_of_dof_per_node; %Total System Degrees of Freedom
assignin('base','Total_System_Degrees_of_Freedom',Total_System_Degrees_of_Freedom);


Degrees_of_Freedom_Per_Element = number_of_nodes_per_element * number_of_dof_per_node; % degrees of freedom per element
assignin('base','Degrees_of_Freedom_Per_Element',Degrees_of_Freedom_Per_Element);

%Fracture Parameters
sigma_t = evalin('base','sigma_t');
Percentage_Limit_Tension = evalin('base','Percentage_Limit_Tension');
tau = evalin('base','tau');
Percentage_Limit_Shear = evalin('base','Percentage_Limit_Shear');
Gf_t = evalin('base','Gf_t');
Percentage_Fracture_Energy_Traction = evalin('base','Percentage_Fracture_Energy_Traction');
d_u = evalin('base','d_u');
Percentage_Ultimate_Deformation_Maximum_Step = evalin('base','Percentage_Ultimate_Deformation_Maximum_Step');
tol = evalin('base','tol');


%This variable will calculate the element size in the
%x-direction
dhx = Length/NXE;

%This variable will calculate the element size in the
%x-direction
dhy = Width/NYE;

X_origin = 0. ; % x origin of the global coordinate system
Y_origin = 0. ; % y origin of the global coordinate system

%Bringing dof from the workspace

is_dx_enabled = evalin('base','is_dx_enabled');
is_dy_enabled = evalin('base','is_dy_enabled');
is_dz_enabled = evalin('base','is_dz_enabled');
is_rx_enabled = evalin('base','is_rx_enabled');
is_ry_enabled = evalin('base','is_ry_enabled');
is_rz_enabled = evalin('base','is_rz_enabled');

%         % This input will ask for the number of degree of freedoms per node
%             %CHANGE number_of_dof_per_element to number_of_dof_per_node

ngpb = evalin('base','ngpb');
ngps = evalin('base','ngps');

%% This Part will ask for NON_GLOBAL Input DATA => That applies to local individual elements

% % This Part will ask for File Inputs necessary for further analysis
% ------------------------------------------------------------------------


% Declaring Important Variables from File
Boundry_Conditions_Degree_of_Freedom = Boundary_Conditions(:,1);%declaring the degree of Freedom of the Boundry Conditions

%%TEMPORARY SPACE Populating the NF matrix-------------
nf = ones(Number_of_Nodes,number_of_dof_per_node); %nodal freedom matrix set to zeros

node_number_where_active = Boundary_Conditions(:,1); %creating a single coumn matrix with node numbers where dof are released
dof_x_displacement = Boundary_Conditions(:,2); %dof release matrix from x-disp column in Boundry condition file
dof_y_displacement = Boundary_Conditions(:,3); %dof release matrix from y-disp column in Boundry condition file
dof_z_displacement = Boundary_Conditions(:,4); %dof release matrix from z-disp column in Boundry condition file
dof_x_rotation = Boundary_Conditions(:,5); %dof release matrix from x-rotation column in Boundry condition file
dof_y_rotation = Boundary_Conditions(:,6); %dof release matrix from y-rotation column in Boundry condition file
dof_z_rotation = Boundary_Conditions(:,7); %dof release matrix from z-rotation column in Boundry condition file
%miss = Boundry_Conditions(:,8); %dof release matrix from last column in Boundry condition file

number_of_active_boundry_conditions = length(node_number_where_active); %total number of nodes at which releases present


%meant to populate the nf matrix
for i = 1:number_of_active_boundry_conditions
    if(dim == 2) %checking if 2D
        if is_dx_enabled
            nf(node_number_where_active(i),1) = dof_x_displacement(i); %Changing nf at x_displacmeent
        end
        if is_dy_enabled
            nf(node_number_where_active(i),2) = dof_y_displacement(i);
        end
        if is_rx_enabled
            nf(node_number_where_active(i),3) = dof_x_rotation(i);
        end
    else
        if(dim ==3) %Checking if 3D
            if is_dx_enabled
                nf(node_number_where_active(i),1) = dof_x_displacement(i); %Changing nf at x_displacmeent
            end
            if is_dy_enabled
                nf(node_number_where_active(i),2) = dof_y_displacement(i); %Changing nf at x_displacmeent
            end
            if is_dz_enabled
                nf(node_number_where_active(i),3) = dof_z_displacement(i); %Changing nf at x_displacmeent
            end
            if is_rx_enabled
                nf(node_number_where_active(i),4) = dof_x_rotation(i); %Changing nf at x_displacmeent
            end
            if is_ry_enabled
                nf(node_number_where_active(i),5) = dof_y_rotation(i); %Changing nf at x_displacmeent
            end
            if is_rz_enabled
                nf(node_number_where_active(i),6) = dof_z_rotation(i); %Changing nf at x_displacmeent
            end
        end
    end
end


assignin('base','nf',nf);

%------------------------------------------------

%meant to count total number of active degrees of
%freedom
total_numbers_of_active_dof = 0;

for i = 1:length(nf)
    for j = 1:number_of_dof_per_node
        if(nf(i,j) == 1) %checking if nf has 1 or 0
            total_numbers_of_active_dof = total_numbers_of_active_dof + 1; %if 1 then increase value by 1
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


%-------------------------------------------------

%for calculation of g in form_KK from platelem_q4

nf_g = nf;
n = 0;
for i=1:Number_of_Nodes
    for j=1:number_of_dof_per_node
        if nf_g(i,j) ~= 0
            n=n+1;
            nf_g(i,j)=n;
        end
    end
end

assignin('base','nf_g',nf_g);

%-------------------------------------------------

%Declaring Important variables from File
% Assign concentrated load


External_load = evalin('base','External_Load');
Nodal_load= External_load(:,1);
Force=External_load(:,2:4);
%
%% Assemble Global force vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1.0 Assign Concentrated load  is the first step and it is ok%
Load = zeros(Number_of_Nodes,3);

for i=1:length(Nodal_load)
    Load(Nodal_load(i),1:3) = Force(i,1:3);
end


assignin('base','Load',Load);
Global_force_vector = zeros(total_numbers_of_active_dof,1);
for i=1:Number_of_Nodes
    for j=1:number_of_dof_per_node
        if nf_g(i,j) ~= 0
            Global_force_vector(nf_g(i,j))= Load(i,j);
        end
    end
end
clear i j
%
% 2.0 Assign gravity load that generates body forces
Gravity_Load=[0 ,-gamma]';

% the gravity load is not related to the gauss point and therefore this part doesn't
% enter in the assign global force vector script

% % 1. Form D matrix of elastic properties
% % [DB]=Dr *[] depends on Dr = Eh3/12(1 − ν2) for bending
% % and  [DS]=G[] depends on G = E/2(1 + ν) for shear
% %
% deeb=formdeeb(Elastic_Modulus,Poissons_Ratio,thickness_of_plate); % Matrix of elastic properties for plate bending
% dees=formdees(Elastic_Modulus,Poissons_Ratio,thickness_of_plate); % Matrix of elastic properties for plate shear
% --------------------------------------------------------------------------
% Input data for Numerical Integration
% --------------------------------------------------------------------------
% % 
% % % 1.Form the matrix containing the abscissas and the weights of Gauss points
% % sampb=gauss(ngpb); samps=gauss(ngps);

% <<<<<<< HEAD
% fg_gravity = zeros(total_numbers_of_active_dof, 1);
% =======
% 1.Form the matrix containing the abscissas and the weights of Gauss points
sampb=gauss(ngpb); samps=gauss(ngps);

% >>>>>>> d15f757c9cc3da402168f8d03084c38dfdd2cd3b

% I discussed with my supervisor regarding this traction load and he
% decides to not approach this part so we can cancel 


current_row = 1;

fg_gravity = fg_matrix_calculator(Number_of_Elements,nodal_connectivity_values,Number_of_Nodes,nf);

Global_force_vector = fg_gravity(1:total_numbers_of_active_dof,1) + Global_force_vector;

assignin('base','fg_gravity',fg_gravity);


%
%% Assemble  global stiffness matrix  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Form D matrix of elastic properties
% [DB]=Dr *[] depends on Dr = Eh3/12(1 − ν2) for bending
% and  [DS]=G[] depends on G = E/2(1 + ν) for shear
%
deeb=formdeeb(Elastic_Modulus,Poissons_Ratio,thickness_of_plate); % Matrix of elastic properties for plate bending
dees=formdees(Elastic_Modulus,Poissons_Ratio,thickness_of_plate); % Matrix of elastic properties for plate shear
%--------------------------------------------------------------------------
% Input data for Numerical Integration
%--------------------------------------------------------------------------

% 1.Form the matrix containing the abscissas and the weights of Gauss points
sampb=gauss(ngpb); samps=gauss(ngps);
% 2. initialize the global stiffness matrix to zero
Global_stiffness_matrix = zeros(total_numbers_of_active_dof,total_numbers_of_active_dof);
%

for iel=1:Number_of_Elements           % loop for the total number of elements
    if Element_Type==3 && ngpb==0
        [bee,fun_3,g,A,d_3] = elem_T3(iel);
        dee=formdeeb(Elastic_Modulus,Poissons_Ratio,thickness_of_plate);
        ke=thickness_of_plate*A*bee'*dee*bee; % Integrate stiffness matrix

        Global_stiffness_matrix=form_KK(Global_stiffness_matrix,ke, g);

    else
        [coord,g] = platelem_q4(iel); % coordinates of the nodes of element i,
        % and its steering vector
        keb=zeros(Degrees_of_Freedom_Per_Element,Degrees_of_Freedom_Per_Element) ; % Initialize the element bending
        % stiffness matrix to zero
        kes=zeros(Degrees_of_Freedom_Per_Element,Degrees_of_Freedom_Per_Element) ; % Initialize the element Shear
        % stiffness matrix to zero

        %--------------------------------------------------------------------------
        % Integrate element Shear stiffness and assemble it in global matrix
        %--------------------------------------------------------------------------
        for ig=1: ngps
            wi = samps(ig,2);
            for jg=1: ngps
                wj=samps(jg,2);
                [der,fun] = fmquad(samps, Element_Type, dim, ig,jg); % Derivative of shape functions
                % in local coordinates
                jac=der'*coord;                    % Compute Jacobian matrix
                d=det(jac);                       % Compute determinant of
                % Jacobian matrix
                jac1=inv(jac);                    % Compute inverse of the
                % Jacobian
                deriv=jac1*der';                   % Derivative of shape functions
                % in global coordinates
                bees=formbees(deriv,fun,number_of_nodes_per_element,Degrees_of_Freedom_Per_Element); % Form matrix [B]
                kes=kes + (5/6)*d*wi*wj*bees'*dees*bees; % Integrate stiffness matrix
            end
        end
        Global_stiffness_matrix=form_KK(Global_stiffness_matrix,kes, g); % assemble global stiffness matrix

    end
    %

end

assignin('base','KK',Global_stiffness_matrix);


%% Displaying Results
%%%%%%%%%%%%%%%%%%%%%%% End of assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp('Solve for unknown displacements' )
%
delta = Global_stiffness_matrix\Global_force_vector;

assignin('base','delta',delta);
%

%

if Element_Type == 3 %Parlas
    for i=1: Number_of_Nodes %
        if nf_g(i,1) == 0 %
            x_disp =0.; %
        else
            x_disp = delta(nf_g(i,1)); %
        end
        %
        if nf_g(i,2) == 0 %
            y_disp = 0.; %
        else
            y_disp = delta(nf_g(i,2)); %
        end
        disp([i x_disp y_disp]) % Display displacements of each node
        DISP(i,:) = [ x_disp y_disp];
    end
end


format short e
disp('node w_disp x_slope y_slope ') %
for i=1: Number_of_Nodes %
    if nf_g(i,1) == 0 %
        w_disp =0.; %
    else
        w_disp = delta(nf_g(i,1)); %
    end
    %
    if nf_g(i,2) == 0 %
        x_slope = 0.; %
    else
        x_slope = delta(nf_g(i,2)); %
    end
    %
    if number_of_dof_per_node > 2
        if nf_g(i,3) == 0 %
            y_slope = 0.; %
            disp([i w_disp x_slope y_slope]) % Display displacements of each node
            DISP(i,:) = [ w_disp x_slope y_slope];
        else
            y_slope = delta(nf_g(i,3)); %
            disp([i w_disp x_slope y_slope]) % Display displacements of each node
            DISP(i,:) = [ w_disp x_slope y_slope];
        end
    else
        disp([i w_disp x_slope]) % Display displacements of each node
        DISP(i,:) = [ w_disp x_slope];
    end


end
%
disp('Calculate moments and shear forces the center of each element')
%
ngp=1; %
samp=gauss(ngp);
%
% Still needs Implementation
for i=1:Number_of_Elements
    [coord,g] = platelem_q4(i); % coordinates of the nodes of element i,
    % and its steering vector
    eld=zeros(Degrees_of_Freedom_Per_Element,1); % Initialize element displacement to zero
    for m=1:Degrees_of_Freedom_Per_Element %
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
            [der,fun] = fmquad(samp, Element_Type, dim, ig,jg,i); % Derivative of shape functions
            % in local coordinates
            jac=der'*coord; % Compute Jacobian matrix
            d=det(jac); % Compute the determinant of
            % Jacobian matrix
            jac1=inv(jac); % Compute inverse of the Jacobian
            deriv=jac1*der'; % Derivative of shape functions
            % in global coordinates
            %
            beeb=formbeeb(deriv,number_of_nodes_per_element,Degrees_of_Freedom_Per_Element); % Form matrix [B_b]
            chi_b = beeb*eld ; % compute bending curvatures
            Moment = deeb*chi_b ; % Compute moments
            bees=formbees(deriv,fun,number_of_nodes_per_element,Degrees_of_Freedom_Per_Element); % Form matrix [B_s]
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


[MX, MY, MXY, QX, QY] = Forces_at_nodes_plate(Element_Forces);

Results = [MX, MY, MXY, QX, QY];
figure;
plot(MX,QX);

assignin('base','Results',Results)

%
figure;
cmin = min(W);
cmax = max(W);
%       caxis([cmin cmax]);
patch('Faces',connec, 'Vertices', geom, 'FaceVertexCData',W,...
    'Facecolor','interp','Marker','.');
colorbar;
toc




