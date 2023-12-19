function [Vector_B, Internal_force_vector] = Vector_B()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%This that guys code


nnd = N(1); %Number of nodes
nel = N(2) ; %Number of elements
nodof = N(3) ; %Number of node dof 
nne = N(4) ; %Number of element nodes
thick = N(5); % Element thickness

% DoF matrix
Dof = zeros(nel,nne*nodof) ;
for i = 1 : nel
    for ii = 1 : nne
        Dof(i,(2*ii-1:2*ii)) = [2*C(i,ii)-1,2*C(i,ii)] ;
    end
end
clear i
clear ii

% Solve the system taking into accoun the BC's
K = zeros(nodof*nnd) ;
f = zeros(nodof*nnd,1) ;
% Assemble the stiffness matrix
 for i = 1 : nel
     [Ke,fe] = stif2d(et,ex(i,:),ey(i,:),thick,b,D) ;
     [K,f] = assem(nne,nodof,C(i,:),K,Ke,f,fe) ;
 end
 clear i
 [u]=solve(K,F,bc) ;
 % Extract element nodal displacements from a global solution vector
 es = []; %element stress matrix
 ep = []; %element strain matrix
 Internal_force_vector = []; %element internal forces matrix
 gpc = []; %gauss points coordinates
 for i = 1 : nel
     ed(i,1:2*nne) = u(Dof(i,:))' ; %transform the nodal disp vector into element disp matrix
     [ess,epp,gpcc]=stress(et,ex(i,:),ey(i,:),D,ed(i,:)') ; %compute stress and strain
     es = [es; ess] ;
     ep = [ep; epp] ;
     gpc = [gpc; gpcc] ;
     if nne==3
         ngp = 1 ;    %triangle number of gauss points        
     elseif nne==4
         ngp = 4 ;    %quad gauss number of points 
     end
     %Computes Internal forces
     [eff]=fint(et,ex(i,:),ey(i,:),thick,es(ngp*(i-1)+1:i*ngp,:)) ;
     Internal_force_vector = [Internal_force_vector; eff] ;
 end
 clear i
       
      
Ft = zeros(size(Cp,1),2);
for i = 1 : size(Cpp,1)
        for ii = Cpp(i,1) : Cpp(i,2)
            for iii = Bee(i,1) : Bee(i,2) 
                for iiii = 1 : nne
                    if Cp(ii) == C(Be(iii),iiii)
                        Ft(ii,:) = [Ft(ii,:) + Internal_force_vector(Be(iii),(2*iiii-1):2*iiii)] ;
                    end
                end
            end
        end
end 

clear i ii iii iiii 


% Vector b  

% Computes the rotation matrix that rotates points in the plane counterclockwise through 
% an angle θ about the local reference system in respect to the global
%two-dimensional Cartesian coordinate system  x1,x2.


if clockwise == 1
    Rotation_matrix = [cos(beta) sin(beta);
            sin(beta) -cos(beta)];
elseif anticlockwise == 2
    Rotation_matrix = [cos(beta) sin(beta);
            -sin(beta) cos(beta)];
end
    

for i = 1 : size(Cp,1)
    RR = Rotation_matrix(2*i-1:2*i,:) ;
    Vector_B(2*i-1:2*i,1) = RR* Ft(i,nodof-1:nodof)' ;
end

    clear i 

end

