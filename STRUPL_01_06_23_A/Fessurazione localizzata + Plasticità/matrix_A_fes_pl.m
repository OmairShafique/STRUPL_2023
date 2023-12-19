function [A, Fnf, Fns] = matrix_A_fes_pl(et,ex,ey,Dof,N,C,bc,exL,eyL,D,Cp,Be,Cpp,Bee,NN,mp,NNN)
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

A = [] ;



for aa = 1 : (size(Cp,1) + size(mp,1))
    
    if aa >= 1 && aa <= size(Cp,1)   
        
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
            
            % Calcolo gli spostamenti dovuti alle distorsioni imposte
            [u]=solve(K,F,bcL) ;
            
            % Calcolo le forze Ft generate dalla distorsione
            Ft = K*u ;
            
            % Calcolo le forze nodali elemento per elemento
            Fnf = zeros(nne*ndof, ne) ; %ogni colonna mi da le forze nodali di ogni elemento
            df = zeros(nne*ndof, ne) ; %ogni colonna mi da gli spostamenti nodali di ogni elemento
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
                        df(:,ww) = d ;
                    end
                end
            end
            
            %Calcolo gli sforzi negli elementi dovuti alle distorsioni imposte
            es_f = []; %element stress matrix
            
            ed_f = df' ;
            
            for j = 1 : ne
                [ess]=stress(et,ex(j,:),ey(j,:),D,ed_f(j,:)') ; %compute stress
                es_f = [es_f; ess] ; %stress
            end
            clear j    
            
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

            % Calcolo le forze Ftt (=-Ft)
            Ftt = K*u ;

            % Calcolo le forze nodali elemento per elemento
            Fns = zeros(nne*ndof, ne) ; %ogni colonna mi da le forze nodali di ogni elemento
            for bb = 1 : size(Cpp,1)            
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
            end

            %Calcolo gli sforzi negli elementi dovuti alle distorsioni imposte
            es_s = []; %element stress matrix
            for j = 1 : ne
                ed_s(j,1:2*nne) = u(Dof(j,:))' ; %transform the nodal disp vector into element disp matrix
            end
            clear j
            for j = 1 : ne
                [ess]=stress(et,ex(j,:),ey(j,:),D,ed_s(j,:)') ; %compute stress
                es_s = [es_s; ess] ; %stress
            end
            clear j
            
            % Somma il contributo a nodi fissi e quello a nodi spostabili
            FTOT = Fnf + Fns  ;
            es_T = es_f + es_s ;
            
            %Proietto le forze
            FNNN = zeros(2*size(Cp,1),1) ;
            for cc = 1 : size(Cp,1)
                RR = NN(2*cc-1:2*cc,:) ; %matrice di rotazione
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
            
            %Proietto gli sforzi
            es_Tp =[] ;
            for j = 1 : size(mp,1)
                i = mp(j) ;
                es_Tpp = es_T(i,:) * NNN(j,:)' ;
                es_Tp = [es_Tp; es_Tpp];
            end
            clear j i
            
            %Assemblo la matrice A
            AA = [-FNNN; -es_Tp];
            A = [A AA] ;
            
        end
        
    else 
        i = aa - size(Cp,1) ;
        j = mp(i) ; %individua a quale elemento è associato il modo plastico
        
        % -1- APPLICO LE DISTORISIONI DIFFUSE NEGLI ELEMENTI A NODI FISSI
        nn_e = NNN(i,:) ; %normale dell'elemento i considerato
        
        %Calcolo lo sforzo nella struttura a nodi fissi dovuto alla
        %deformazione nn_n
        es_ff = -D*nn_e' ;
        es_f = es_ff' ;
        
        %Matrice degli sforzi della struttura
        es_F = zeros(ne,3);
        es_F(j,:) = es_f ;
        
        %Le forze ef sono le forze nodali nell'elemento i dovute allo
        %sforzo es_f
        [ef] = fint(et,ex(j,:),ey(j,:),t,es_f) ;
        
        %Considerando la struttura a nodi fissi, l'unico elemento con sforzi
        %diversi da zero (e quindi anche forze nodali diverse da zero) è quello
        %in cui impongo le deformazioni.
        
        %Assemblo una matrice con le forze nodali 
        Fnf = zeros(nne*ndof, ne) ; %ogni colonna mi da le forze nodali di ogni elemento
        Fnf(:,j) = ef' ;
        
        %Calcolo quindi le forze nodali Fint (assemblate in un vettore globale),
        %cambio segno e applico alla struttura a nodi spostabili.
        
        Fint = zeros(ndof*nn,1) ;
        
        for jj = 1 : size(C,2)
            k = C(j,jj) ;
            Fint(2*k-1) = ef(2*jj-1);
            Fint(2*k) = ef(2*jj);
        end
        clear jj
             
        % -2- RISOLVO LA STRUTTURA A NODI SPOSTABILI
        [u_n]=solve(K,-Fint,bc) ;
        
        % Calcolo le forze nodali elemento per elemento
        Fns = zeros(nne*ndof, ne) ; %ogni colonna mi da le forze nodali di ogni elemento
        for bb = 1 : size(Cpp,1)
            for dd = Bee(bb,1) : Bee(bb,2)
                CC = C(Be(dd),:) ; %nodi dell'elemento considerato
                [Ke] = stif2d(et,exL(dd,:),eyL(dd,:),t,b,D) ;
                d = zeros(ndof*nne,1) ;
                for ii = 1 : nne %Number of element nodes
                    for jj = 1 : ndof %Number of node dof
                        rw = ndof*(CC(ii)-1)+jj ; %Row position
                        d(ndof*(ii-1)+jj)= u_n(rw); %Displacement of element nodes
                    end
                end
                clear ii jj
                ff = Ke*d ;  %Forze sui nodi dell'elemento considerato
                ww = Be(dd) ;
                Fns(:,ww) = ff ;
            end
        end

        %Calcolo gli sforzi a partire dagli spostamenti
        es_s = []; %element stress matrix
        for jj = 1 : ne
            ed_n(jj,1:2*nne) = u_n(Dof(jj,:))' ; %transform the nodal disp vector into element disp matrix
        end
        clear jj
        for jj = 1 : ne
            [ess]=stress(et,ex(jj,:),ey(jj,:),D,ed_n(jj,:)') ; %compute stress
            es_s = [es_s; ess] ; %stress
        end
        clear jj

        % Somma il contributo a nodi fissi e quello a nodi spostabili
        FTOT = Fnf + Fns  ;
        es_T = es_F + es_s ;
        
        %Proietto le forze
        FNNN = zeros(2*size(Cp,1),1) ;
        for cc = 1 : size(Cp,1)
            RR = NN(2*cc-1:2*cc,:) ; %matrice di rotazione
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
        
        %Proietto gli sforzi
        es_Tp =[] ;
        for jj = 1 : size(mp,1)
            ii = mp(jj) ;
            es_Tpp = es_T(ii,:) * NNN(jj,:)' ;
            es_Tp = [es_Tp; es_Tpp];
        end
        clear jj ii
        
        %Assemblo la matrice A
        AA = [-FNNN; -es_Tp];
        A = [A AA] ;
    
 end

end