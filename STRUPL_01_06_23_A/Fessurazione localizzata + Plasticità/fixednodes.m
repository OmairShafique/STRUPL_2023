function [FT, FNF,Df] = fixednodes(Lambda,exL,eyL,C,Cp,Cpp,Be,Bee,N,NN,et,D)
%FIXEDNODES
% La funzione risolve la struttura a nodi fissi soggetta alle distorsioni
% imposte (lambda)

% Structure data
nn = N(1); %Number of nodes
ne = N(2) ;
ndof = N(3);  %Number of node dof 
nne = N(4) ; %Number of element nodes
t = N(5); % Element thickness

% The structures is considered unloaded
F = zeros(nn*ndof,1); %Global load vector
b = [0  0]' ; %Element body force vector 

FT = zeros(nn*ndof,1);
FNF = zeros(nne*ndof, ne) ;
Df = zeros(nne*ndof, ne) ;

 for aa = 1 : size(Lambda,1)  
    % -1- RISOLUZIONE DELLA STRUTTURA A NODI FISSI CON  DISTORSIONE UNITARIA
    % Applico una distorsione unitaria in direzione n sui nodi del crack path
    nx = NN(aa,1) ;
    ny = NN(aa,2) ;
    
    if rem(aa,2)==0 
        aaa = aa/2 ;
        r = Cp(aaa,1) ;
    else
        aaa = (aa+1)/2 ;
        r = Cp(aaa,1) ;
    end
    
    bcL = [1:1:(ndof*nn);zeros(1,ndof*nn)]' ;
    bcL(ndof*r-1,2) = -Lambda(aa)*nx ;
    bcL(ndof*r,2) = -Lambda(aa)*ny ;

        % Solving the system taking into accoun the BC's
        K = zeros(ndof*nn) ;
        f = zeros(ndof*nn,1) ;
        % Assembling of the lower part stiffness matrix
        for bb = 1 : size(Cpp,1)
            if (aaa >= Cpp(bb,1)) && (aaa <= Cpp(bb,2))
                for i = Bee(bb,1) : Bee(bb,2) 
                    [Ke,fe] = stif2d(et,exL(i,:),eyL(i,:),t,b,D) ;
                    [K,f] = assem(nne,ndof,C(Be(i),:),K,Ke,f,fe) ;
                end
            end
        end
        clear i bb
        
        % Calcolo le forze Ft generate ai nodi dalla distorsione
        [u]=solve(K,F,bcL) ;
        % Calcolo le forze generate ai nodi elemento per elemento
        Fnf = zeros(nne*ndof, ne) ; %ogni colonna mi da le forze nodali di ogni elemento
        df = zeros(nne*ndof, ne) ;
        for bb = 1 : size(Cpp,1)
             if (aaa >= Cpp(bb,1)) && (aaa <= Cpp(bb,2)) %verifico a quale linea di fe
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
                        ff = Ke*d ;  %Forze sui nodi dell'elemento considerato 
                        ww = Be(dd) ;
                        Fnf(:,ww) = ff ;
                        df(:,ww) = d ;
                    end 
             end
        end
        
        Df = Df + df ;

        %Calcolo le forze nodali sommando il contributo dei diversi
        %elementi
        Ft  = zeros(ndof*nn,1) ;
        for bb = 1 : size(Cpp,1)
            if (aaa >= Cpp(bb,1)) && (aaa <= Cpp(bb,2)) %verifico a quale linea di fessura appartiene il nodo aa
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

        %Calcolo delle forze totali
        FT = FT + Ft ;
        FNF = FNF + Fnf ;
 end
end

