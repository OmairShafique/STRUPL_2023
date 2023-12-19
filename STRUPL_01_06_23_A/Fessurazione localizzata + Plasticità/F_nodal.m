function [ F_n] = F_nodal( Cpp, Bee, Be, C, et, exL, eyL, t, b, D, ndof, nne, ne, u )
%F_nodal 
%Calcola, a partire dagli spostamenti nodali della struttura u, le forze 
%nodali F_n elemento per elemento.
%Solo le colonne relative agli elementi di influenza dei nodi di fessura
%avranno valori diversi da zero
%-------------------------------------------------------------
% INPUT:    
%           u    :  Nodal displacement (nn*ndof,1)
%           C    :  Element matrix
%           ndof :  Number of dof per node
%           nne  :  Number of nodes per element
%           nn   :  Number of nodes of the structure
%           ne   :  Number of elements of the structure
%           Cpp  :  Distribution of cracking nodes in cracking lines
%           Be   :  Boundouray elements 
%           Bee  :  Distribution of boundary elements on cracking lines
% 
% OUTPUT:   F_projection    :  Force orojection on the normal and
%                              tangential directions (nne*ndof,ne)

%--------------------------------------------------------------
    F_n = zeros(nne*ndof,ne) ;
    for aa = 1 : size(Cpp,1)
        for bb = Bee(aa,1) : Bee(aa,2) 
            CC = C(Be(bb),:) ; %nodi dell'elemento considerato
            [Ke] = stif2d(et,exL(bb,:),eyL(bb,:),t,b,D) ;
            d = zeros(ndof*nne,1) ;
            for ii = 1 : nne %Number of element nodes
                for jj = 1 : ndof %Number of node dof
                    rw = ndof*(CC(ii)-1)+jj ; %Row position     
                    d(ndof*(ii-1)+jj)= u(rw); %Displacement of element nodes
                end
            end
            clear ii jj
            ff = Ke*d ;  %Forze sui nodi dell'elemento considerato 
            cc = Be(bb) ;
            F_n(:,cc) = ff ;
        end 
    end

end

