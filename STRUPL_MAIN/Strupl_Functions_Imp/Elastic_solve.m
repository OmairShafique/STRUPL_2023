function[delta,Reaction,node_displacement,SIGMA,STRAIN]=Elastic_solve()


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
for i=1:nnd;
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