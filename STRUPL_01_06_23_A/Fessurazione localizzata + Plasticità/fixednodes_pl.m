function [es] = fixednodes_pl(Fint,K,Dof,ex,ey,N,et,D)
%FIXEDNODES
% La funzione risolve la struttura a nodi fissi soggetta alle distorsioni
% imposte (lambda)

% Structure data
nn = N(1); %Number of nodes
ne = N(2) ;
ndof = N(3);  %Number of node dof 
nne = N(4) ; %Number of element nodes
t = N(5); % Element thickness


b = [0  0]' ; %Element body force vector 
bc = [1:1:(ndof*nn);zeros(1,ndof*nn)]' ; %Boundary conditions

% Calcolo gli spostamenti generati dalla distorsione
[u]=solve(K,Fint,bc) ;

% Extract element nodal displacements from a global solution vector
es = []; %element stress matrix
ep = []; %element strain matrix

for i = 1 : ne
    ed(i,1:2*nne) = u(Dof(i,:))' ; %transform the nodal disp vector into element disp matrix  
end
clear i

for i = 1 : ne
    [ess,epp]=stress(et,ex(i,:),ey(i,:),D,ed(i,:)') ; %compute stress and strain
    es = [es; ess] ; %stress
    ep = [ep; epp] ; %strain
end
clear i

 end

