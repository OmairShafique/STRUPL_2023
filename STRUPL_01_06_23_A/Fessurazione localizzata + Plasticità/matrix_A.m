function [A, Fnf, Fns] = matrix_A(et,ex,ey,N,C,bc,exL,eyL,D,Cp,Be,Cpp,Bee,NN)
%matrix_A calculates the A matrix in order to solve the PLCP.
%Both w (normal crack) and s (tangetial crack) are taking into account 
%in order to be able to cath the collapse mechanism when a piecewise 
%potential crack is considered
%-------------------------------------------------------------
% INPUT:    ex   :  Nodes coordinates in x direction
%           ey   :  Nodes coordinates  in y direction
%           et   :  Element type: 1=4 node quad 2=3 node triangle
%           C    :  Element matrix
%           b    :  Body force vector in x and y direction 
%           D    :  Constitutive matrix
%           u    :  Nodal displacements
%           N    =  [nn ne ndof nne t]
%           nn   :  Total number of nodes
%           ne   :  Total numer of elements
%           ndof :  Number of dof per node
%           nne  :  Number of nodes per element
%           ex   :  Nodes coordinates  in x direction
%           ey   :  Nodes coordinates  in y direction
%           Cp   :  Potential cracking nodes
%           Cpp  :  Distribution of cracking nodes in cracking lines
%           Be   :  Boundouray elements for each cracking line
%           Bee  :  Distribution of boundary elements on cracking lines
%           BETA :  Angles between x axes and the normal to the crack lines
% 
% OUTPUT:   A    :  Matrix A
%           FT   :  Reazioni a nodi fissi
%--------------------------------------------------------------

nn = N(1); %Number of nodes
ne = N(2) ; %Number of elements
ndof = N(3) ; %Number of node dof 
nne = N(4) ; %Number of element nodes
t = N(5); % Element thickness

% The structures is considered unloaded
F = zeros(nn*ndof,1); %Global load vector
b = [0  0]' ; %Element body force vector 

 for aa = 1 : size(Cp,1)  

    % -1- RISOLUZIONE DELLA STRUTTURA A NODI FISSI CON  DISTORSIONE UNITARIA
    % Applico una distorsione unitaria in direzione n sui nodi del crack path
    nx = NN(2*aa-1,1) ;
    ny = NN(2*aa-1,2) ;
    tx = NN(2*aa,1) ;
    ty = NN(2*aa,2) ;
      
    for pp = 1 : 2 %applico la distorsione prima in direzione normale e poi tangenziale
        %A seconda della linea a cui appartiene il nodo selezionato, individuo
        %il rispettivo angolo beta per calcolare le componenti della distorsione
        if pp == 1     
        % Applico una distorsione unitaria in direzione normale alla fessura
            r = Cp(aa,1) ;
            bcL = [1:1:(ndof*nn);zeros(1,ndof*nn)]' ;
            bcL(ndof*r-1,2) = -nx;
            bcL(ndof*r,2) = -ny ;
        elseif pp == 2
            % Applico una distorsione unitaria in direzione tangenziale alla fessura
            r = Cp(aa,1) ;
            bcL = [1:1:(ndof*nn);zeros(1,ndof*nn)]' ;
            bcL(ndof*r-1,2) =  -tx;
            bcL(ndof*r,2) = -ty ;
        end
        
        % Solving the system taking into accoun the BC's
        K = zeros(ndof*nn) ;
        f = zeros(ndof*nn,1) ;
        % Assembling of the lower part stiffness matrix
        for bb = 1 : size(Cpp,1)
            if (aa >= Cpp(bb,1)) && (aa <= Cpp(bb,2))
                for i = Bee(bb,1) : Bee(bb,2) 
                    [Ke,fe] = stif2d(et,exL(i,:),eyL(i,:),t,b,D) ;
                    [K,f] = assem(nne,ndof,C(Be(i),:),K,Ke,f,fe) ;
                end
            end
        end
        clear i bb
        
        % Calcolo le forze Ft generate ai nodi dalla distorsione
        [u]=solve(K,F,bcL) ;

        Ft = K*u ;
        % Calcolo le forze generate ai nodi elemento per elemento
        Fnf = zeros(nne*ndof, ne) ; %ogni colonna mi da le forze nodali di ogni elemento
        for bb = 1 : size(Cpp,1)
            if (aa >= Cpp(bb,1)) && (aa <= Cpp(bb,2)) %verifico a quale linea di fessura appartiene il nodo aa
                    for dd = Bee(bb,1) : Bee(bb,2) 
                        CC = C(Be(dd),:) ; %nodi dell'elemento considerato
                        [Ke] = stif2d(et,exL(dd,:),eyL(dd,:),t,b,D) ;
                        d = zeros(ndof*nne,1) ;
                        for ii = 1 : nne %Number of element nodes
                            for jj = 1 : ndof %Number of node dof
                                rw = ndof*(CC(ii)-1)+jj ; %Row position     
                                d(ndof*(ii-1)+jj)= u(rw) ;%Displacement of element nodes
                            end
                        end
                        clear ii jj
                        ff = Ke*d  ; %Forze sui nodi dell'elemento considerato 
                        ww = Be(dd) ;
                        Fnf(:,ww) = ff ;
                    end 
            end
        end

        %Calcolo le forze nodali sommando il contributo dei diversi
        %elementi
        Ft  = zeros(ndof*nn,1) ;
        for bb = 1 : size(Cpp,1)
            if (aa >= Cpp(bb,1)) && (aa <= Cpp(bb,2)) %verifico a quale linea di fessura appartiene il nodo aa
                    for dd = Bee(bb,1) : Bee(bb,2) 
                        for ii = 1 : nne
                            for jj = 1 : ndof 
                                CC = C(Be(dd),ii) ;  %numero del nodo considerato nell'elemento dd
                                ww = Be(dd) ;
                                Ft(ndof*CC+(jj-2),1) = Ft(ndof*CC+(jj-2),1)+ Fnf(ndof*ii+(jj-2),ww) ;
                            end
                        end
                    end
            end
        end
    
% -2- RISOLUZIONE DELLA STRUTTURA GLOBALE SOGGETTA A -Ft  (NODI SPOSTABILI)
        F = Ft ;
        % Assembling of the global stiffness matrix
        K = zeros(ndof*nn) ;
        f = zeros(ndof*nn,1) ;
        for i = 1 : ne
            [Ke,fe] = stif2d(et,ex(i,:),ey(i,:),t,b,D) ;
            [K,f] = assem(nne,ndof,C(i,:),K,Ke,f,fe) ;
        end
        clear i
    
        % Calcolo gli spostamenti della struttura a cui ho applicato -Ft
        [u]=solve(K,-F,bc) ;
        Ftt = K*u ;
        % Calcolo le forze generate ai nodi elemento per elemento
        Fns = zeros(nne*ndof, ne) ; %ogni colonna mi da le forze nodali di ogni elemento
        for bb = 1 : size(Cpp,1)
%              if (aa >= Cpp(bb,1)) && (aa <= Cpp(bb,2)) %verifico a quale linea di fessura appartiene il nodo aa
                    for dd = Bee(bb,1) : Bee(bb,2) 
                        CC = C(Be(dd),:) ; %nodi dell'elemento considerato
                        [Ke] = stif2d(et,exL(dd,:),eyL(dd,:),t,b,D) ;
                        d = zeros(ndof*nne,1) ;
                        for ii = 1 : nne %Number of element nodes
                            for jj = 1 : ndof %Number of node dof
                                rw = ndof*(CC(ii)-1)+jj ; %Row position     
                                d(ndof*(ii-1)+jj)= u(rw); %Displacement of element nodes
                            end
                        end
                        clear ii jj
                        ff = Ke*d ;  %Forze sui nodi dell'elemento considerato 
                        ww = Be(dd) ;
                        Fns(:,ww) = ff ;
                    end 
%             end
        end
        
        % Somma il contributo a nodi fissi e quello a nodi spostabili
        FTOT = Fnf + Fns  ;
        
        RR = [nx ny; tx ty] ; %matrice di rotazione
        
        FNNN = zeros(2*size(Cp,1),1) ;
        for cc = 1 : size(Cp,1) 
            RR = NN(2*cc-1:2*cc,:) ;
            for bb = 1 : size(Cpp,1)
                if (cc >= Cpp(bb,1)) && (cc <= Cpp(bb,2)) %verifico a quale linea di fessura appartiene il nodo cc
                    for dd = Bee(bb,1) : Bee(bb,2) 
                        for ii = 1 : nne
                            CC = C(Be(dd),ii) ;  %numero del nodo considerato nell'elemento dd
                            ww = Be(dd) ;
                            if Cp(cc) == CC
                                FNNN(2*cc-1:2*cc)= FNNN(2*cc-1:2*cc) + RR * FTOT(ndof*ii-1:ndof*ii,ww) ;
                            end
                        end
                    end
                end
            end
        end                     
        if aa==1 && pp==1
            A = [-FNNN] ;
        else
            A = [A -FNNN] ;    
        end
               
    end
    
 end

end