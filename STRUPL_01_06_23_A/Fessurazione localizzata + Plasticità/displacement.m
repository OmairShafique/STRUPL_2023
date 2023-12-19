function [ed,u,u_t] = displacement(alpha,F,K,bc,Lambda,exL,eyL,C,Cp,Cpp,Be,Bee,N,NN,et,D,Dof)
%Calcola gli spostamenti della struttura elemento per elemento%
%-------------------------------------------------------------
% INPUT:    
%           alpha    :  Load multiplier
%           F        :  External loads
%           K        :  Stiffness matrix
%           bc       :  Boundary conditions
%           Lambda   :  Plastic deformations in local reference system
%           exL      :  Boundary Element X coordinates
%           eyL      :  Boundary Element Y coordinates
%           C        :  Element matrix
%           Cp       :  Potential cracking nodes
%           Cpp      :  Distribution of cracking nodes in cracking lines
%           Be       :  Boundouray elements 
%           Bee      :  Distribution of boundary elements on cracking lines
%           N        :  Structur data N = [nn ne ndof nne t]
%           NN       :  Rotation matrix (size(Cp,1)*ndof,ndof)
%           et       :  Element type: 1=4 node quad 2=3 node triangle
%           D        :  Elastic matrix 
%           Dof      :  DoF matrix
% 
% OUTPUT:   u        :  Compatible displacement vector (nn*ndof,1)
%           u_t      :  Plastic deformations in global reference system (nn*ndof,1)
%           ed       : Total displacement(element by element)(nne*ndof,ne)

%--------------------------------------------------------------
%Calcolo la soluzione elastica della struttura soggetta a Fc
        nn = N(1); %Number of nodes
        ne = N(2) ; %Number of elements
        ndof = N(3); %Number of node dof 
        nne = N(4); %Number of element nodes
        Fc = alpha*F  ;
        [u_c]=solve(K,Fc,bc);  %
        %Calcolo le forze Ft generate da Lambda sulla struttura a nodi fissi
        [Ft] = fixednodes(Lambda,exL,eyL,C,Cp,Cpp,Be,Bee,N,NN,et,D) ;
        %Calcolo la soluzione elastica della struttura soggetta a -Ft
        [u_p]=solve(K,-Ft,bc) ; 
        %Calcolo lo spostamento totale 
        u = u_c + u_p;  %compatible displacements
        u_t = zeros(nn*ndof,1) ;    
        %Proiezione delle lambda su x e y
        for ccc = 1 : size(Cp,1) 

            RRR = NN(2*ccc-1:2*ccc,:) ;
            RR = RRR(:,1)' ;
            ux(ccc,1) = RR * Lambda(2*ccc-1:2*ccc) ;
            RR = RRR(:,2)' ;
            uy(ccc,1) = RR * Lambda(2*ccc-1:2*ccc) ;
        end

        u_t(2*Cp-1) = ux ;    
        u_t(2*Cp) = uy ; 

%        Extract element nodal displacements from a global solution vector
 es = []; %element stress matrix
 ep = []; %element strain matrix

  for i = 1 : ne
     ed(i,1:2*nne) = u(Dof(i,:))' ; %transform the nodal disp vector into element disp matrix  
  end
  d = zeros(size(ed));
  for cc = 1 : size(Cp,1)
      for bb = 1 : size(Cpp,1)
          if (cc >= Cpp(bb,1)) && (cc <= Cpp(bb,2))
              for dd = Bee(bb,1) : Bee(bb,2)
                  for ii = 1 : nne
                      CC = C(Be(dd),ii) ;   % nodo elemento considerato 
                      ww = Be(dd) ;
                      if Cp(cc) == CC
                          d(ww,2*ii-1) = d(ww,2*ii-1)  - ux(cc) ;
                          d(ww,2*ii) = d(ww,2*ii) - uy(cc)  ;
                      end
                  end
              end
          end
      end
  end
  ed = ed + d ;
                  
end

