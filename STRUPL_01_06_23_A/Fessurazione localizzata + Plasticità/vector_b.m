function [bb, ef] = vector_b(et,ex,ey,N,C,F,b,bc,D,Cp,Cpp,Be,Bee,NN)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nn = N(1); %Number of nodes
ne = N(2) ; %Number of elements
ndof = N(3) ; %Number of node dof 
nne = N(4) ; %Number of element nodes
t = N(5); % Element thickness

% DoF matrix
Dof = zeros(ne,nne*ndof) ;
for i = 1 : ne
    for ii = 1 : nne
        Dof(i,(2*ii-1:2*ii)) = [2*C(i,ii)-1,2*C(i,ii)] ;
    end
end
clear i
clear ii

% Solve the system taking into accoun the BC's
K = zeros(ndof*nn) ;
f = zeros(ndof*nn,1) ;
% Assemble the stiffness matrix
 for i = 1 : ne
     [Ke,fe] = stif2d(et,ex(i,:),ey(i,:),t,b,D) ;
     [K,f] = assem(nne,ndof,C(i,:),K,Ke,f,fe) ;
 end
 clear i
 [u]=solve(K,F,bc) ;
 % Extract element nodal displacements from a global solution vector
 es = []; %element stress matrix
 ep = []; %element strain matrix
 ef = []; %element internal forces matrix
 gpc = []; %gauss points coordinates
 for i = 1 : ne
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
     [eff]=fint(et,ex(i,:),ey(i,:),t,es(ngp*(i-1)+1:i*ngp,:)) ;
     ef = [ef; eff] ;
 end
 clear i
       
      
Ft = zeros(size(Cp,1),2);
for i = 1 : size(Cpp,1)
        for ii = Cpp(i,1) : Cpp(i,2)
            for iii = Bee(i,1) : Bee(i,2) 
                for iiii = 1 : nne
                    if Cp(ii) == C(Be(iii),iiii)
                        Ft(ii,:) = [Ft(ii,:) + ef(Be(iii),(2*iiii-1):2*iiii)] ;
                    end
                end
            end
        end
end 

clear i ii iii iiii 

% Vector b      
for i = 1 : size(Cp,1)
    RR = NN(2*i-1:2*i,:) ;
    bb(2*i-1:2*i,1) = RR* Ft(i,ndof-1:ndof)' ;
end

    clear i 

end

