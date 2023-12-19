function [ F_p ] = F_projection( F, ndof, nne, NN, Cp, Cpp, Be, Bee, C )
%F_projection 
%Calcola, a partire dalle forze nodali F di ogni elemento, la
%proiezione normale e tangenziale nei nodi di potenziale fessura Cp,
%tenendo conto degli opportuni elementi di influenza Be e delle matrice
%di rotazione.
%-------------------------------------------------------------
% INPUT:    
%           F    :  Nodal force matrix (nne*ndof,size(Be,1))
%           C    :  Element matrix
%           ndof :  Number of dof per node
%           nne  :  Number of nodes per element
%           ne   :  Number of elements of the structure
%           Cp   :  Potential cracking nodes
%           Cpp  :  Distribution of cracking nodes in cracking lines
%           Be   :  Boundouray elements 
%           Bee  :  Distribution of boundary elements on cracking lines
%           NN   :  Rotation matrix (size(Cp,1)*ndof,ndof)
% 
% OUTPUT:   F_projection    :  Force orojection on the normal and
%                              tangential directions

%--------------------------------------------------------------
F_p = zeros(2*size(Cp,1),1) ;
    for aa = 1 : size(Cp,1)
        RR = NN(2*aa-1:2*aa,:) ;
        for bb = 1 : size(Cpp,1)
            if (aa >= Cpp(bb,1)) && (aa <= Cpp(bb,2)) %verifico a quale linea di fessura appartiene il nodo aa
                for cc = Bee(bb,1) : Bee(bb,2) 
                    for dd = 1 : nne
                        CC = C(Be(cc),dd) ;  %numero del nodo considerato nell'elemento dd
                        ww = Be(cc) ;
                        if Cp(aa) == CC
                            F_p(2*aa-1:2*aa)= F_p(2*aa-1:2*aa) + RR * F(ndof*dd-1:ndof*dd,ww) ;
                        end
                    end
                end
            end
        end
    end

end

