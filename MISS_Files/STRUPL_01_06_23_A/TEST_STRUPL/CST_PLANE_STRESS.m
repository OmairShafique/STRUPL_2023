% THIS PROGRAM USES AN 3-NODE LINEAR TRIANGULAR ELEMENT FOR THE
% LINEAR ELASTIC STATIC ANALYSIS OF A TWO DIMENSIONAL PROBLEM
%
clear all
clc
tic
%
% Make these variables global so they can be shared by other functions
%
global nnd nel nne nodof eldof n
global geom connec dee nf Nodal_loads
%
format long g
%
% ALTER NEXT LINES TO CHOOSE THE NAME OF THE OUTPUT FILE
%
fid =fopen('CST_COARSE_MESH_RESULTS.txt','w');
%
% To change the size of the problem or change elastic properties
% supply another input file
%
CST_COARSE_MESH_DATA;
%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Assemble the global force vector
%
fg=zeros(n,1);
for i=1: nnd
if nf(i,1) ~= 0
fg(nf(i,1))= Nodal_loads(i,1);
end
if nf(i,2) ~= 0
fg(nf(i,2))= Nodal_loads(i,2);
end
end
%
% Assembly of the global stiffness matrix
%
% initialize the global stiffness matrix to zero
%
KK = zeros(n, n);
%
for i=1:nel
[bee,g,A] = elem_T3(i); % Form strain matrix, and steering vector
ke=thick*A*bee'*dee*bee; % Compute stiffness matrix
KK=form_KK(KK,ke, g); % assemble global stiffness matrix
end
%
%
%%%%%%%%%%%%%%%%%%%%%%% End of assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
delta = KK\fg ; % solve for unknown displacements
%
node_disp=zeros(nnd,2);
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
node_disp(i,:) =[x_disp y_disp];
end
%
% Retrieve the x_coord and y_disp of the nodes located on the neutral axis
%
k = 0;
for i=1:nnd;
if geom(i,2)== 0.
k=k+1;
x_coord(k) = geom(i,1);
vertical_disp(k)=node_disp(i,2);
end
end
%
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
EPS(i,:)=eps ; % Store strains for all elements
sigma=dee*eps; % Compute stresses
SIGMA(i,:)=sigma ; % Store strains for all elements
end
%
% Print results to file
%
print_CST_results;
%
% Plot the stresses in the x_direction
%
x_stress = SIGMA(:,1);
cmin = min(x_stress);
cmax = max(x_stress);
clim([cmin cmax])
patch('Faces', connec, 'Vertices', geom, 'FaceVertexCData',x_stress, ...
'Facecolor','flat','Marker','o')
colorbar;

toc