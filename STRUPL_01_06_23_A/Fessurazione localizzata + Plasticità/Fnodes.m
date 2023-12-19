function [Fr] = Fnodes(exL,eyL,et,Be,Bee,N,C,Cp,Cpp,b,D,u)
%Fnodes calculates the nodal forces Fr due to the nodal compatible 
%displacment u of some elements (BE) of the structure 
%-------------------------------------------------------------
% INPUT:    ex   :  Nodes coordinates in x direction
%           ey   :  Nodes coordinates  in y direction
%           et   : Element type: 1=4 node quad 2=3 node triangle
%           BE   :  Elements considered (total or part of the structure=
%           C    :  Element matrix
%           b    :  Body force vector in x and y direction 
%           D    :  Constitutive matrix
%           u    :  Nodal displacements
%           N    =  [nn ne ndof nne t]
%           nn   :  Total number of nodes
%           ne   :  Total numer of elements
%           ndof :  Number of dof per node
%           nne  :  Number of nodes per element
% 
% OUTPUT:   Fr   : Nodal forces 
%--------------------------------------------------------------

nn = N(1); %Number of nodes
ne = N(2) ; %Number of elements
ndof = N(3) ; %Number of node dof 
nne = N(4) ; %Number of element nodes
t = N(5); % Element thickness

Fr = zeros(ndof*nn,1) ;
    for bb = 1 : size(Cpp,1) %linea di fessura
        for cc = Cpp(bb,1) : Cpp(bb,2) %nodi sulla linea di fessura
            for dd = Bee(bb,1) : Bee(bb,2) %elementi linea di fessura
                CC = C(Be(dd),:) ;
                for ee = 1 : size(CC,2) %nodi dell'elemento
                    if Cp(cc) == CC(ee) %
                        [Ke,fe] = stif2d(et,exL(dd,:),eyL(dd,:),t,b,D) ;
                        d = zeros(ndof*nne,1) ;
                        for ii = 1 : nne %Number of element nodes
                            for jj = 1 : ndof %Number of node dof
                                rw = ndof*(CC(ii)-1)+jj ; %Row position        
                                d(ndof*(ii-1)+jj)= u(rw) ;%Displacement of element nodes
                            end
                        end
                        clear ii jj
                        ff = Ke*d ;  %Forze sui nodi dell'elemento considerato                  
                        for ii = 1 : nne
                            for jj = 1: ndof
                                if  CC(ii) == CC(ee) 
                                    %Considero solo le forze relative al nodo considerato
                                    Fr(ndof*(CC(ii)-1)+jj) = Fr(ndof*(CC(ii)-1)+jj)+ ...
                                    ff(ndof*(ii-1)+jj) ;
                                end
                            end
                        end
                        clear ii jj d ff Ke fe 
                    end
                end
            end
        end   
    end
end

