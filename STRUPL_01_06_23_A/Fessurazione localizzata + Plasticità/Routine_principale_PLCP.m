clear all
close all
clc
cputime

% Risoluzione problema di collasso strutturale tramite PLCP
% Elementi finiti numerati sempre in verso antiorario
% Per la risoluzione numerica si utilizzano le trasformazioni pivotali
% ATTENZIONE! NUMERO MAX MODI ATTIVABILI DEVE ESSERE < 10000"
% Elementi da considerare rispetto alla linea di fessura: inferiori


% !!!!!!!!!!!!!!!!!!!! ATTENZIONE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% SE VENGONO VIOLATI I VINCOLI DOPO SCARICO CAMBIARE MIN CON MAX
% O VICEVERSA (RIGA 980)

%-------------------------------------------------------------------------
%------------------------------ INPUT DATA -------------------------------
%-------------------------------------------------------------------------

tol = 1.e-4 ;

% Definizione del problema
% 0 = Olonomo
% 1 = Anolonomo
TP = 1 ;

%Number of step
step = 50;

%Geometry
load('Coordinates.txt');
load('Element_matrix.txt');
Cd = Coordinates ;
C = Element_matrix ;
ngp = 1;

%Elements coordinates
for i = 1 : size(C,1)
    ex(i,:) = Cd(C(i,:),1);  %Element X coordinates
    ey(i,:) = Cd(C(i,:),2);  %Element Y coordinates
end
clear i

%Structure data
nn = size(Cd,1); %Number of nodes
ndof = size(Cd,2) ; %Number of node dof
ne = size(C,1) ; %Number of elements
nne = size(C,2) ; %Number of element nodes
t = 420 ; % Element thickness
N = [nn ne ndof nne t];

%Material data
at = 1 ; % Analysis type: 1=plane stress  2=plane strain
et = 2 ; % Element type: 1=4 node quad 2=3 node triangle
nu = 0.2 ;  % Poisson ratio
E = 1500; % Elastic modulus MPa

%Macro-cracks
sigma_t = 0.1198; % Limit tension [MPa]
tau = 0.1198; % Limit shear [MPa]
tau_m = tau*3 ; %Valore di tau per macrofess
Gf_t = 0.035; % Fracture Energy traction [N/mm]
Gf_s = Gf_t ; % Fracture Energy shear[N/mm]
w_t = 2 * Gf_t/sigma_t ;
w_s = 2 * Gf_s/(tau) ;
h_s_t = -sigma_t /w_t  ; %softening parameter
h_s_s = -tau /w_s  ; %softening parameter

%Micro-cracks (Drucker-Prager parameters)
mp = [1:1:ne]' ; %modi plastici
d_u = 0.019 ; %Ultimate deformation
coeff_h = 0.01 ; %percentuale di resistanza in corrispondenza del quale si attiva l'hardning
% h_t = ((1-coeff_h)*sigma_t /d_u) ; %hardening parameter
% h_s = ((1-coeff_h)*tau/d_u) ; %hardening parameter
h_t = (0.9*sigma_t /d_u) ; %hardening parameter


sigma_t_r = coeff_h * sigma_t ; % Limit traction
sigma_c = sigma_t_r  ; % Limit compression
%-------------------------
sigma_t_tol = sigma_t_r + eps*sigma_t_r; 
sigma_c_tol = sigma_c + eps*sigma_c ;
Aa = (2/sqrt(3))*((sigma_t_r*sigma_c)/(sigma_t_r+sigma_c)) ;
Bb = (1/sqrt(3))*((sigma_t_r-sigma_c)/(sigma_t_r+sigma_c)) ;
Aa_tol = (2/sqrt(3))*((sigma_t_tol*sigma_c_tol)/(sigma_t_tol+sigma_c_tol)) ;
Bb_tol = (1/sqrt(3))*((sigma_t_tol-sigma_c_tol)/(sigma_t_tol+sigma_c_tol)) ;



%Volume of elements
Vol = [];
for i = 1 : mp
    Vol = [Vol; polyarea(ex(mp(i),:),ey(mp(i),:))*t] ;
end

% DoF matrix
Dof = zeros(ne,nne*ndof) ;
for i = 1 : ne
    for ii = 1 : nne
        Dof(i,(2*ii-1:2*ii)) = [2*C(i,ii)-1,2*C(i,ii)] ;
    end
end
clear i
clear ii

%Loads and BC's
b = [0  0]' ; % Element body force vector ef
load('External_load.txt');
F = External_load(:,2) ;
% Fx = zeros(nn,1);
% Fy = zeros(nn,1);
% Fxy = External_load ;
% for ff = 1 : size(Fxy,1)
%     Fx(Fxy(ff,1)) = Fxy(ff,2);
%     Fy(Fxy(ff,1)) = Fxy(ff,3);
% end
% F = zeros(2*nn,1);


load('Boundary_condition.txt');
bc = Boundary_condition ;

% Elastic matrix
D = hooke(at,E,nu) ;

%Stiffness matrix
K = zeros(ndof*nn) ;
f = zeros(ndof*nn,1) ;
for i = 1 : ne
    [Ke,fe] = stif2d(et,ex(i,:),ey(i,:),t,b,D) ;
    [K,f] = assem(nne,ndof,C(i,:),K,Ke,f,fe) ;
end
clear i

%Cracks data
load('Beta.txt');
load('Cracking_node.txt');
Cp = Cracking_node' ;
load('Cracking_node_distribution.txt');
Cpp = Cracking_node_distribution ;
load('Boundary_element.txt');
Be = Boundary_element' ;
load('Boundary_element_distribution.txt');
Bee = Boundary_element_distribution ;

%Boundary elements coordinates
for i = 1 : size(Be,1)
    exL(i,:) = Cd(C(Be(i),:),1);  %Element X coordinates
    eyL(i,:) = Cd(C(Be(i),:),2);  %Element Y coordinates
end

%Assemblage of nodal rotation matrix
for cc = 1 : size(Cp,1)
    if cc == 1
        [NN] = Rotation_matrix(Beta(cc,2),Beta(cc,3)) ;
    else
        [NN] = [NN; Rotation_matrix(Beta(cc,2),Beta(cc,3))] ;
    end
end


%--------------------------------------------------------------------------
%Calcolo delle normali per i modi plastici
%--------------------------------------------------------------------------
% Risolvo il problema elastico per carico assegnato
[ui]=solve(K,F,bc) ;

% Extract element nodal displacements from a global solution vector
es_o = []; %element stress matrix
for i = 1 : ne
    %transform the nodal disp vector into element disp matrix
    ed(i,1:2*nne) = ui(Dof(i,:))' ;
end
clear i j
for i = 1 : ne
    [ess]=stress(et,ex(i,:),ey(i,:),D,ed(i,:)') ; %compute stress
    es_o = [es_o; ess] ; %stress
end
clear i j

%Calcolo per ogni elemento il moltiplicatore che attiva il modo
for j = 1 : size(mp,1)
    i = mp(j) ;
    %Calcolo delle tensioni principali
    es_po(i,1) = ((es_o(i,1)+es_o(i,2))/2) + (((es_o(i,1)-es_o(i,2))/2)^2 + (es_o(i,3)^2))^(1/2) ; %sigmaI
    es_po(i,2) = ((es_o(i,1)+es_o(i,2))/2) - (((es_o(i,1)-es_o(i,2))/2)^2 + (es_o(i,3)^2))^(1/2) ; %sigmaII
    es_po(i,3) = (((es_o(i,1)-es_o(i,2))/2)^2 + (es_o(i,3)^2))^(1/2) ; %taumax
    
    %Calcolo invariante primo
    I1(j,1) = es_po(i,1) + es_po(i,2) ;
    
    %Calcolo invariante deviatorico secondo
    J2(j,1) = 1/6*(es_o(i,1)-es_o(i,2))^2 + 1/6*es_o(i,1)^2 + 1/6*es_o(i,2)^2 + es_o(i,3)^2 ;
    
    alpha_b(j,1) = Aa/(sqrt(J2(j))-Bb*I1(j)) ;
end
clear i j

%Calcolo le normali
ed = []; %element strain matrix
NNN = [] ;
ES = [] ;
for j = 1 : size(mp,1)
    % Risolvo il problema elastico per carico amplificato di alpha
    [u_alpha]=solve(K,alpha_b(j)*F,bc) ;
    % Extract element nodal displacements from a global solution vector
    es = []; %element stress matrix
    
    for i = 1 : ne
        %transform the nodal disp vector into element disp matrix
        ed(i,1:2*nne) = u_alpha(Dof(i,:))' ;
    end
    clear i
    
    for i = 1 : ne
        [ess]=stress(et,ex(i,:),ey(i,:),D,ed(i,:)') ; %compute stress
        es = [es; ess] ; %stress
    end
    clear i
    
    i = mp(j) ;
    ES = [ES; es(i,:)] ;
    
    %Calcolo invariante deviatorico secondo
    J2_ = 1/6*(es(i,1)-es(i,2))^2 + 1/6*es(i,1)^2 + 1/6*es(i,2)^2 + es(i,3)^2 ;
    
    %Calcolo della normale in un punto della superficie limite
    dphi_dsigmax = -Bb + 0.5*J2_^(-0.5) * (1/3*(es(i,1)-es(i,2))+1/3*es(i,1)) ;
    dphi_dsigmay = -Bb + 0.5*J2_^(-0.5) * (-1/3*(es(i,1)-es(i,2))+1/3*es(i,2)) ;
    dphi_dsigmaxy =  0.5*J2_^(-0.5) * 2*es(i,3) ;
    
    nnn = [dphi_dsigmax dphi_dsigmay dphi_dsigmaxy]  ;
    nnn_n = nnn/norm(nnn) ;
    NNN = [NNN; nnn_n] ;
end
clear j



%--------------------------------------------------------------------------
%---------------------PLCP Vectors and matrix------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% VETTORE DELLE RESISTENZE ------------------------------------------------
%--------------------------------------------------------------------------
% Il vettore delle resistenze è costituito da R_f che contiene le
% resistenze relative ai modi di fessurazione localizzata (prodotto della
% tensione  limite per l'area di influenza di ciascun nodo, calcolata a
% fuori dal codice) e R_p che contiene la resistenza relativa ai modi
% plastici (prodotto dello sforzo dell'elemento per la rispettiva normale)
% La resistenza associata ai modi di softening in R_f è zero.


load('Influence_lenght.txt'); %lunghezze di influenza ai nodi
IL = Influence_lenght(:,2) ;

% ILL = Influence_lenght(:,2) ;
% IL = zeros(2*size(ILL,1),2) ;
% IL(1:2:end)= ILL(1:end) ;
% IL(2:2:end)= ILL(1:end) ;

RR = zeros(size(IL)) ;
RR(1:2:end) = IL(1:2:end)*t*sigma_t;
RR(2:2:end) = IL(2:2:end)*t*tau_m;
 
% RR(108) = RR(108)*1.5;

% RR(117) = RR(117)*0.5;
% 
% RR(113) = RR(113)*1.5;
% RR(121) = RR(121)*1.5;
% 
% RR(173) = RR(173)*10;
% RR(247) = RR(247)*10;
% RR(227) = RR(227)*10;
% RR(193) = RR(193)*10;9
% RR(17) = RR(17)*2;
% RR(61) = RR(61)*2;


R_f = zeros(2*2*size(Cp,1),1) ;
for rr = 1 : size(RR,1)
    R_f(2*rr-1) = RR(rr) ;
end
clear rr

R_p = [] ;
for j = 1 : size(mp,1)
    i = mp(j) ;
    RR_p = ES(j,:)* NNN(j,:)' ;
    R_p = [R_p; RR_p] ;
end
clear j
%Moltiplico per il volume di ciscun elemento (intragrazione nel volume)
R_p = R_p.*Vol ;

R = [R_f; R_p];


%--------------------------------------------------------------------------
% VETTORE DEI CARICHI b ---------------------------------------------------
%--------------------------------------------------------------------------
% La componente relativa ai modi plastici bb_p si calcola moltiplicando gli
% sforzi elastici dovuti ai carichi esterni per la rispettiva normale

bb_p = [] ;
for j = 1 : size(mp,1)
    es_op = es_o(mp(j),:) * NNN(mp(j),:)' ;
    bb_p = [bb_p; es_op] ;
end
clear j
%Moltiplico per il volume di ciscun elemento (intragrazione nel volume)
bb_p = bb_p.*Vol ; 

% La componente associata alla fessurazione localizzata bb_f si calcola con
% apposita funzione e viene poi modificata per tenere conto del softening

[bbf, ef] = vector_b(et,ex,ey,N,C,F,b,bc,D,Cp,Cpp,Be,Bee,NN) ;

bb_f = zeros(2*size(bbf,1),1) ;
for xx = 1 : size(bbf,1)
    bb_f(2*xx-1) = bbf(xx);
    bb_f(2*xx) = bbf(xx);
end
clear xx

bb = [bb_f; bb_p] ;


%--------------------------------------------------------------------------
% MATRICE A ---------------------------------------------------------------
%--------------------------------------------------------------------------

[AA] = matrix_A_fes_pl(et,ex,ey,Dof,N,C,bc,exL,eyL,D,Cp,Be,Cpp,Bee,NN,mp,NNN) ;

% Scompongo la matrice A in quadranti e li modifico opportunamente per
% tenere conto dell'hardening dei modi plastici e del softening dei modi
% della fessura localizzata

% Quadrante plastico
for yy = 2*size(Cp,1) + 1 :  size(AA,1)
    for zz = 2*size(Cp,1) + 1 :  size(AA,1)
        yyy = yy - 2*size(Cp,1) ;
        zzz = zz - 2*size(Cp,1) ;
        if yyy == zzz 
            A_p(yyy,zzz) = AA(yy,zz) + h_t ;
        else
            A_p(yyy,zzz) = AA(yy,zz) ;
        end
    end
end
clear yy yyy zz zzz
for i = 1 : size(A_p,2) 
    %Moltiplico per il volume di ciscun elemento (intragrazione nel volume)
    A_p(:,i) = A_p(:,i).*Vol ;
end
clear i

% Quadrante fessura localizzata
H = zeros(2*size(Cp,1),1);
for hh = 1 : size(H,1)
    if rem(hh,2)==1 
        H(hh,1) = h_s_t * t * IL(hh);
    else
        H(hh,1) = h_s_s * t * IL(hh);
    end
end
clear hh

for yy = 1 : 2*size(Cp,1)
    for zz = 1 : 2*size(Cp,1)
         if yy == zz
            A_f(2*yy-1,2*zz-1) = AA(yy,zz)+ H(yy);
            A_f(2*yy,2*zz) = AA(yy,zz);
            A_f(2*yy-1,2*zz) = AA(yy,zz) + H(yy);
            A_f(2*yy,2*zz-1) = AA(yy,zz);
%         elseif rem(yy,2)==1 && zz == yy + 1
%             A_f(2*yy-1,2*zz-1) = AA(yy,zz)+ H(yy);
%             A_f(2*yy,2*zz) = AA(yy,zz);
%             A_f(2*yy-1,2*zz) = AA(yy,zz) + H(yy);
%             A_f(2*yy,2*zz-1) = AA(yy,zz);  
%         elseif rem(yy,2)==0 && zz == yy - 1
%             A_f(2*yy-1,2*zz-1) = AA(yy,zz)+ H(yy);
%             A_f(2*yy,2*zz) = AA(yy,zz);
%             A_f(2*yy-1,2*zz) = AA(yy,zz) + H(yy);
%             A_f(2*yy,2*zz-1) = AA(yy,zz); 
        else
            A_f(2*yy-1,2*zz-1) = AA(yy,zz);
            A_f(2*yy,2*zz) = AA(yy,zz);
            A_f(2*yy-1,2*zz) = AA(yy,zz);
            A_f(2*yy,2*zz-1) = AA(yy,zz);
        end
    end
end
clear yy zz

A_f_rid = A_f(1:2:end, 1:2:end);

% Quadrante fessura localizzata + plastico
for yy = 2*size(Cp,1) + 1 :  size(AA,1)
    for zz = 1 : 2*size(Cp,1)
        xx = yy - 2*size(Cp,1) ;
        A_fp(xx,2*zz-1) = AA(yy,zz) ;
        A_fp(xx,2*zz) = AA(yy,zz) ;
    end
end
clear xx yy zz

for i = 1 : size(A_fp,2)
    %Moltiplico per il volume di ciscun elemento (intragrazione nel volume)
    A_fp(:,i) = A_fp(:,i).*Vol ; 
end
clear i

A_fp_rid = A_fp(:,1:2:end);

% Quadrante plastico + fessura localizzata
for yy = 1 : 2*size(Cp,1)
    for zz = 2*size(Cp,1) + 1 :  size(AA,1)
        xx = zz - 2*size(Cp,1) ;
        A_pf(2*yy-1,xx) = AA(yy,zz) ;
        A_pf(2*yy,xx) = AA(yy,zz) ;
    end
end
clear xx yy zz

A_pf_rid = A_pf(1:2:end,:);

A = [A_f A_pf; A_fp A_p] ;

A_rid = [A_f_rid A_pf_rid; A_fp_rid A_p] ;


%--------------------------------------------------------------------------
%------------------------ Plastic problem solution ------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% STEP 0 - SOLUZIONE AL LIMITE ELASTICO  ----------------------------------
%--------------------------------------------------------------------------

% Costruzione del Tableau iniziale T --------------------------------------

%Softening
P_f = zeros(size(A_f,1),1) ;
P_f(1:2:end) = 10000*(1:1:2*size(Cp,1)) ;
P_f(2:2:end) = -10000*(1:1:2*size(Cp,1)) ;
L_f = (P_f/10000)';

P_p = zeros(size(A_p,1),1) ;
P_p(1:1:end) = 10000*(2*size(Cp,1)+1 : 1 :(size(mp,1)+2*size(Cp,1))) ;
L_p = (P_p/10000)';

L = [L_f L_p];
P = [P_f; P_p];

T = [0  100000000  200000000 L ; P R -bb A] ;
Tinitial = T ;
To = T ;
To(:,2:3) = [] ; % La matrice A con le intestazioni (senza colonna delle
%resisteze e dei carichi) mi serve per l'unloading

% % Aggiungo la riga e la colonna di intestazione(base e fuori base)per
% % facilitare gli scambi del Tableau utilizzando dei numeri identificativi
% % scelti arbitrariamente.
% % La corrispondenza tra lambda e phi è rappresentata da un fattore 10000.
% % lambda = 1,-1,2,-2...
% % phi = 10000,-10000,200,-200...
% % R = 100000000
% % alpha = 200000000
% % I modi col meno sono quelli relativi al softening


% Ricerca dell'elemento Pivot ---------------------------------------------

Alpha = T(:,2)./T(:,3) ;
p = max(Alpha(Alpha<0)) ;
cc = 3 ; %colonna elemento Pivot (nota in partenza)
rr = find(Alpha==p) ; %riga elemento Pivot
if size(rr,1)>1
    rr = rr(1);
end
rro = rr;
cco = cc;

% Costruzione del nuovo tableau T_new -------------------------------------

% Il nuovo tableau si costruisce applicando le regole di pivoting
% Il primo elemento ad entrare in base è sempre Alpha
% Si utilizza la funzione 'pivoting' che effettua la trasformazione
% restituendo il nuovo tableauàà

[T_new] = pivoting( T,rr,cc );
Telastic = T_new;
alpha0 = T_new((find((T_new(:,1))==200000000)),2) ;  %200000000=identifica alpha
act_mod = T_new(1,cc) ;

% Costruzione dei vettori Phi e Lambda ------------------------------------

Phi0 = zeros(size(P)) ;
for pp = 1 : size(P,1)
    for tt = 1 : size(P,1)
        if P(pp,1)== T_new(tt+1,1)
            Phi0(pp)=T_new(tt+1,2) ;
        end
    end
end
clear ff
clear tt

Lambda0 = zeros(size(L')) ;
for ll = 1 : size(L,2)
    for tt = 1 : size(L,2)
        if L(1,ll)== T_new(tt+1,1)
            Lambda0(ll)=T_new(tt+1,2) ;
        end
    end
end

% Cacolo spostamenti e sforzi ---------------------------------------------

Fa = alpha0*F ; % forze al limite elastico
[u_el, R_el]=solve(K,Fa,bc) ; % spostamenti al limite eleastico
[Fcc] = Fnodes(exL,eyL,et,Be,Bee,N,C,Cp,Cpp,b,D,u_el) ;

%Calcolo le forze nodali negli elementi di influenza a partire dagli
%spostamenti
[F_n] = F_nodal(Cpp,Bee,Be,C,et,exL,eyL,t,b,D,ndof,nne,ne,u_el) ;
%Proietto le forze tenendo conto degli elementi di influenza al nodo
[Fne] = F_projection(F_n,ndof,nne,NN,Cp,Cpp,Be,Bee,C) ;
clear ccc

% Calcolo gli sforzi e gli invarianti
es_el = []; %element stress matrix
ep_el = []; %element strain matrix
for i = 1 : ne
    ed(i,1:2*nne) = u_el(Dof(i,:))' ; %transform the nodal disp vector into element disp matrix
    [ess,epp]=stress(et,ex(i,:),ey(i,:),D,ed(i,:)') ; %compute stress and strain
    es_el = [es_el; ess] ;
    ep_el = [ep_el; epp] ;
    %Calcolo invariante primo
    I1_el(i,1) = es_el(i,1) + es_el(i,2) ;
    %Calcolo invariante deviatorico secondo
    J2_el(i,1) = 1/6*(es_el(i,1)-es_el(i,2))^2 + 1/6*es_el(i,1)^2 + 1/6*es_el(i,2)^2 + es_el(i,3)^2 ;
end
clear ess epp

% Controllo che i punti siano sulla superficie di plasticità (proietto
% sforzi)
es_elP = [] ;
    for i = 1 : size(mp,1)
        j = mp(i) ;
        es_elp = es_el(j,:) * NNN(i,:)';
        es_elP = [es_elP; es_elp];
    end
    clear j i  
    
U_TOT = [u_el] ;
RF = [R_el] ;
FN = [Fne ] ;
ES_TOT = [es_el] ;
ES_TOTP =[es_elP] ;
J2_TOT = [J2_el] ;

%-------------------------------------------------------------------------- 
% Inizializzaione vettori e check
alpha_s = 0 ; %diventa diverso da zero quando si attiva uno scarico
Lambda_s = 0 ; %diventa diverso da zero quando si attiva uno scarico
LAMBDAS = [];
ET = zeros(step,1) ; %vettore delle trasformazioni vuote
Tcc = 0 ; %vettore per memorizzare i modi eliminati dal tableau
Tcheck = 0 ; %check sul tipo di trasformazione
Tp = 0 ;
checkload = 0 ; %check sull'andamento del carico (crescende/decrescente)
Acheck = 0 ; %check sulla matrice A (def. pos/negativa)
Un_olo = 0 ;
tncheck = 0 ; %check sul termino noto della matrice alla Lemke
unload = 0 ; %check su scarichi avvenuti
inv = 0; %check sull'inversione di carico
T_La = 0 ;
collapse = 0 ;
vincoli = 0 ;
 
%-------------------------------------------------------------------------- 
% ------------------------ Pivoting successivi ----------------------------
%-------------------------------------------------------------------------- 

for ii = 1 : step % numero step da scegliere di volta in volta
    S = [' step ', num2str(ii)];
    disp(S)
    
    T = T_new ;
    %----------------------------------------------------------------------
    % Individuo la nuova colonna e riga pivot -----------------------------
    %----------------------------------------------------------------------

    if T(1,cc)<10000 && T(1,cc)>0 && T(rr,1)==T(1,cc)*10000 
        %C'è stata una trasfrmazione vuota allo step precedente
        display 'LambdaB entra in base'
        Tcheck(ii) = 1 ;
        if Tcc == 0
            Tcc = T(1,cc) ;
        else
            Tcc = [Tcc; T(1,cc)] ;
        end
        
        i = find(prev_mod==T(rr,1)) ;
        prev_mod(i) = [] ;
        
        T(:,cc)=[] ;
        T(rr,:)=[] ;
        TTTT = T ;      
        cc = find(T(1,:)==-Tcc(end)) ;
        mm= T(:,2)./T(:,cc) ;
        mm(1)= 100000000 ;  %riga di intestazione (metto numero positivo)
        p = max(mm(mm<0)) ;
        rr = find(mm==p);  %riga elemento Pivot
        if size(rr,1)>1
            rr = rr(1);
        end

    elseif T(1,cc) <= -10000 %l'ultimo modo attivato è PhiB(<10000)
        display 'Trasformazione vuota'
        Tcheck(ii) = 2 ;
        %eseguo la trasformazione vuota scambiando PhiA e LambdaA
        cc = find(T(1,:)==(-T(1,cc))) 
        rr = find(T(:,1)==(T(1,cc)/10000)) 
        if size(rr,1)>1
            rr = rr(1);
        end
        Lambda_b = T(rr,2) ; %memorizzo il valore di LamdaA
        T(rr,2) = 0 ; %azzero il valore di LambdaA
        ET(ii) = 1;
        TTT = T ;
    elseif  T(1,cc) >= 10000 % l'ultimo modo attivato è PhiA(>10000)
        display 'LambdaA entra in base'
        Tcheck(ii) = 3 ;
        TTT = T ;
        cc = find(T(1,:)==T(1,cc)/10000) ;% scambio LambdaA (carico crescente)
        cc
        mm= T(:,2)./T(:,cc) ;
        mm(1)= 100000 ;  %riga di intestazione (metto numero positivo)
        
        %Inserisco un ciclo di controllo nel caso in cui avessi uno zero
        %nella colonna dei termini noti con un termine negativo nella
        %colonna pivot. Senza il check il rapporto uguale a zero non viene
        %considerato.
        for j =  1 : size(mm, 1)
            if T(j,2) == 0 && T(j,cc)<0
                rr = j ;
                break
                
            else
                p = max(mm(mm<0)) ;
                rr = find(mm==p) ; %riga elemento Pivot
            end
        end
        
        
        rr
        %Se dovessi avere due modi con lo stesso valore scelgo uno dei due
        if size(rr,1)>1
            rr = rr(1);
        end
        
        if T(rr,1)< 0 && T(rr,1)>-10000
            display 'LambdaB esce dalla base - cambio riga'
            %Se un modo LambdaB vuole uscire dalla base,seleziono una diversa riga pivot.
            check_LambdaB = 0 ;  %il check mi serve per interrompere il ciclo quando necessario
            while T(rr,1)< 0 && T(rr,1)>-10000
                check_LambdaB = 1 ;
                mm(rr)=0 ; %passo a selezionare una nuova riga
                p = max(mm(mm<0)) ;
                if isempty(p) == 1
                    display 'Errore rapporto > 0 LambdaB - unbounded solution'
                    Err = 1 ;
                    break
                end
                rr = find(mm==p) ;  %riga elemento Pivot
                if size(rr,1)>1
                    rr = rr(1);
                end
                if check_LambdaB == 0
                    break
                else
                    check_LambdaB = 0; %resetto il check
                end
            end
        end

        if T(rr,1) >= 10000
            display 'PhiA esce dalla base'
        elseif T(rr,1) > 0 && T(rr,1) < 10000
            display 'LambdaA esce dalla base'
        elseif T(rr,1) <= -10000
            display 'PhiB esce dalla base'
        end     
        
        rr
        
    else % l'ultimo modo attivato è Lambda(A o B)
        display 'Phi(A o B) entra in base'
        Tcheck(ii) = 4 ;
        cc = find(T(1,:)==T(1,cc)*10000) ; %scambio Lambda e  rispettivo Phi
        display 'scarico'
        mm= T(:,2)./T(:,cc) ;
        mm(1)= 0 ;  %riga di intestazione (metto numero positivo)
        p = max(mm(mm<0)) ;
        rr = find(mm==p);  %riga elemento Pivot
        if size(rr,1)>1
            rr = rr(1);
        end
    end
    
    if T(2:end,cc)>=-tol 
        if T(find(T(:,1)==200000000),cc)<tol && T(find(T(:,1)==200000000),cc)>-tol
            display 'collasso o moltiplicatore zero'
        
%         %Calcolo gli spostamenti di meccanismo
%         Lambda_mec = zeros(size(L')) ;
%         for ll = 1 : size(L,2)
%             for tt = 1 : size(T,1)-1
%                 if L(1,ll)== T(1,cc)
%                     Lambda_mec(ll)=1 ;
%                 elseif L(1,ll)== T(tt+1,1)
%                     Lambda_mec(ll)=T(tt+1,cc) ;
%                 end
%             end
%         end
%         kk = (size(Lambda_mec,1) - size(mp,1))/2 ;
%         LambdaS_mec = zeros(kk + size(mp,1),1) ;
%         for jj = 1 : kk
%             if Lambda_mec(2*jj) > tol
%                 LambdaS_mec(jj) = Lambda_mec(2*jj-1) + Lambda_mec(2*jj) ;
%             else
%                 LambdaS_mec(jj) = Lambda_mec(2*jj-1) ;
%             end
%         end
%         for jj = 2*kk+1 : size(Lambda_mec,1)
%             LambdaS_mec(jj-kk) = Lambda_mec(jj) ;
%         end
%         LambdaS_mecp = LambdaS(1:kk) ;
%         [ed,u,u_t] = displacement(alpha,F,K,bc,LambdaS_mecp,exL,eyL,C,Cp,Cpp,Be,Bee,N,NN,et,D,Dof) ;
        
        break
        end
    end
    
    %----------------------------------------------------------------------
    % Controlli sulla riga pivot
    %----------------------------------------------------------------------
    %%%%%%%I DUE CONTROLLI DOVREBBERO ESSERE CONTESTUALI!!!!!!%%%%%%%%%%%%%
    
    %Se un modo PhiB vuole uscire dalla base, devo prima verificare che
    %il relativo modo PhiA sia fuori dalla base. Se così non fosse,
    %passo a selezionare una diversa riga pivot.
    check = 0 ;  %il check mi serve per interrompere il ciclo quando necessario
  
    while T(rr,1)<= -10000 || T(rr,1) < 0 && T(rr,1)>-10000 || T(rr,1) == 200000000 
        if T(rr,1)<= -10000
            for ff = 1 : size(T,1)
                if T(rr,1)== -T(ff,1)
                    check = 1 ; %il relativo modo PhiA è ancora in base
                    mm(rr)=0 ; %passo a selezionare una nuova riga
                elseif T(rr,1)== -T(1,cc)
                    check = 1 ; %il relativo modo PhiA vuole entrare in base
                    mm(rr)=0 ; %passo a selezionare una nuova riga
                end
            end
        elseif T(rr,1) < 0 && T(rr,1)>-10000 || T(rr,1) == 200000000  %LambdaB e alpha
            check = 1 ;
            mm(rr) = 0 ;
        end
        
        p = max(mm(mm<0)) ;
        if isempty(p) == 1
            display 'errore rapporto>0 - unbounded solution'
            break
        end
        rr = find(mm==p) ;
        %riga elemento Pivot
        if size(rr,1)>1
            rr = rr(1);
        end
        
        if check == 0
            break
        else
            check = 0; %resetto il check
        end
    end
       
    if T(rr,1) >= 10000 && Tcheck(ii) ~= 2
        display 'PhiA esce dalla base'
    elseif T(rr,1) > 0 && T(rr,1) < 10000  && Tcheck(ii) ~= 2
        display 'LambdaA esce dalla base'
    elseif T(rr,1) <= -10000 && Tcheck(ii) ~= 2
        display 'PhiB esce dalla base'
    end
    
    %----------------------------------------------------------------------
    %Memorizzo in dei vettori le posizioni degli elementi pivot
    if ii == 1
        RRRR = [rro rr] ;
        CCC = [cco cc] ;
    elseif ii>1
        RRRR = [RRRR rr] ;
        CCC = [CCC cc] ;
    end
    %----------------------------------------------------------------------
    % Controllo se il problema è semi definito positivo
    %----------------------------------------------------------------------
    if ii > 1
        if Tcheck(ii-1) == 2 %nel caso ci sia stata una trasformazione vuota
            DAlpha = ALPHA(ii) - ALPHA(ii-2) ;
        else
            DAlpha = ALPHA(ii) - ALPHA(ii-1) ;
        end
        DAlpha
        inv
        
        %Nel caso in cui c'è stata una trasformazione vuota, non aggiorno
        %la T_La in modo da avere contemporanemte nel tableau sia il modo
        %olonomo che quello anolonomo
        if Tcheck(ii) ~= 1 
            [T_La,Ap] = ACheck_pl_fes( prev_mod, last_mod, L_p, A,bb);
        end

        [pos, AAA,v,AAA_v,AAA_b] = def_pos_test( prev_mod, last_mod, A, L_p, bb) ;
        if v > 0
            display 'v > 0'
        end
        
        if pos == 2
            display 'Matrice non definita '
        elseif pos == 1
            display 'Determinante nullo - Meccanismo '
        elseif pos == 0
            display 'Matrice definita positiva '
        end
        coeff = T_La(end,end) 
        
        if Tcheck(ii-1) ~= 2
            if DAlpha > 0
                inv = 0 ;
            elseif DAlpha < 0
                inv = 1 ;
            end
        end
        
       
        if TP == 1 %problema anolonomo
            %             if  T(rr_alpha,cc) < 0 && T(rr_alpha,2) > 0
            for ww = 1 : size(T,1)
                if T(ww,1)<= L_p(end) && T(ww,1)>0 && Tcheck(ii)~= 2
                    %Individuo le righe relative ai LambdaA
                    %se il relativo coefficiente della colonna pivot è
                    %negativo significa che avrò uno scarico e la matrice
                    %non è definita positiva, se è positivo vuol dire che
                    %non avrò nessuno scarico locale
                    if T(ww,cc) < 0
                        Acheck = 1 ;
                    end
                end
            end
            %             end
        end
        
        TTTT = T;
    end
    
    %----------------------------------------------------------------------
    % Caso di probelma non semidefinito positivo
    %----------------------------------------------------------------------
    if Acheck == 1     
         %------------------ INIZIO FUNZIONE UNLOADING ---------------------
        %La funzione può essere raccolta in uno script a parte per snellire il
        %codice. Per un maggior controllo della procedura è stata lasciata
        %nel codice.
        
        prev_mod = sort(prev_mod) ;
        Err = 0 ;
        T_La_L = T_La ;

        if last_mod <0
            %Se l'ultimo modo attivato è di tipo olonomo
            [T_A] = Act_mod( prev_mod, last_mod, A, L_p, bb ) ;
            [T_A_nh] = Act_mod_no_olo(T_A) ;
            [Tz_finn,Err, Tz_in] = LCP_lemke_tot(T_A_nh) ;
            Err
            if Err == 1
                T_AA = T_A_nh ;
                %Se non ho trovato soluzione per carico crescente
                %inverto il carico e risolvo LCP per carico decrescente
                T_AA(2:end,2) = -T_AA(2:end,2) ;
                display 'Inversione di carico - LCP'
                inv = 1 ;
                [Tz_finn,Err] = LCP_lemke_tot(T_AA) ;
                Err
            end
            T_La_Lnew = Tz_finn ;
            
            %--------------------------------------------------------------
            unload_mod = [] ; %vettore dei modi che si scaricano
            act_mod = [] ; %vettore dei modi che rimangono attivi
            
            % Costruisco il vettore dei modi che scaricano e uno di quelli che
            % rimangono attivi
            for vv = 3 : size(T_La_Lnew,2)
                %si pescano dalla prima riga i modi ancora attivi;
                if T_La_Lnew(1,vv)>=10000 && T_La_Lnew(1,vv)<50000000  || T_La_Lnew(1,vv)<= - 10000
                    act_mod = [act_mod; T_La_Lnew(1,vv)] ;
                end
            end
            clear vv
            for vv = 2 : size(T_La_Lnew,1) %si pescano dalla prima colonna i modi che scaricano;
                if T_La_Lnew(vv,1)>=10000
                    unload_mod = [unload_mod; T_La_Lnew(vv,1)];
                end
            end
            clear vv
            
            mod = [act_mod; unload_mod] ;
            act_mod
            unload_mod =  sort(unload_mod)
            
        else
            
            %%%%% CARICO CRESCENTE ----------------------------------------
            if inv == 0
                display 'Carico crescente - PLCP'
                [T_La_Lnew, Err] = PLCP_lastmod(T_La_L, last_mod) ;
                T_La_Lnew1 = T_La_Lnew ;
                Err
                if Err == 1
                    %Se non ho trovato soluzione al PLCP per carico crescente
                    %inverto il carico e risolvo LCP per carico decrescente
                    T_La_L(2:end,2) = -T_La_L(2:end,2) ;
                    display 'Inversione di carico - LCP'
                    inv = 1 ;
                    [Tz_fin,Err] = LCP_lemke(T_La_L) ;
                    Err
                    display 'Carico decrescente - PLCP'
                    [T_La_Lnew, Err, T_newp] = PLCP_lastmod(Tz_fin, last_mod) ;
                end
                %%%%% CARICO DECRESCENTE ----------------------------------
            else
                %Se sono in carico decrescente inverto il segno del termine
                %noto e risolvo PLCP
                T_La_L(2:end,2) = -T_La_L(2:end,2) ;
                display 'Carico decrescente - PLCP'
                [T_La_Lnew, Err] = PLCP_lastmod(T_La_L, last_mod) ;
                T_La_Lnew1 = T_La_Lnew ;
                Err
                if Err == 1
                    %Se non ho trovato soluzione PLCP per carico decrescente
                    %inverto il carico e risolvo LCP per carico crescente
                    T_La_L(2:end,2) = -T_La_L(2:end,2) ;
                    display 'Inversione di carico - LCP'
                    inv = 0 ;
                    Err = 0 ;
                    [Tz_fin,Err,Tz_1] = LCP_lemke(T_La_L) ;
                    display 'Carico crescente - PLCP'
                    [T_La_Lnew, Err] = PLCP_lastmod(Tz_fin, last_mod) ;
                end
            end
            
            inv
            
            %%%%%%-----------------------------------------------------
            
            unload_mod = [] ; %vettore dei modi che si scaricano
            act_mod = [] ; %vettore dei modi che rimangono attivi
            
            % Costruisco il vettore dei modi che scaricano e uno di quelli che
            % rimangono attivi
            for vv = 3 : size(T_La_Lnew,2)
                %si pescano dalla prima riga i modi ancora attivi;
                if T_La_Lnew(1,vv)>=10000 || T_La_Lnew(1,vv)<= - 10000
                    act_mod = [act_mod; T_La_Lnew(1,vv)] ;
                end
            end
            clear vv
            for vv = 2 : size(T_La_Lnew,1) %si pescano dalla prima colonna i modi che scaricano;
                if T_La_Lnew(vv,1)>=10000
                    unload_mod = [unload_mod; T_La_Lnew(vv,1)];
                end
            end
            clear vv
            
            mod = [act_mod; unload_mod] ;
            
            act_mod
            unload_mod = sort(unload_mod)
        end
        
        %%%%%%-------------------------------------------------------------
        
        %Richiamo il tableau iniziale per ripartire e sostituisco il valore di zero
        %nella colonna dei termini noti all'alpha e a tutte le lambda che erano attive prima
        %dello scarico
        T_new = T ;
        for jj = 1 : size(T,1)
            if T(jj,1)>0 && T(jj,1)<10000
                T_new(jj,2) = 0 ;
            elseif T(jj,1)<0 && T(jj,1)>-10000
                T_new(jj,2) = 0 ;
            elseif T(jj,1) == 200000000
                T_new(jj,2) = 0 ;
            end
        end
        T0=T_new ;
        unload = 1 ;

       
        %Effettuo le trasformazioni vuote per portare in base le Phi
        %che si sono scaricate
        if unload_mod ~= 0
            for zz = 1 :  size(unload_mod,1)
                T = T_new ;
                cc = find(T(1,:)==unload_mod(zz,1)) ;  %colonna elemento Pivot
                rr = find(T(:,1)==unload_mod(zz,1)/10000) ; %riga elemento Pivot
                [T_new] = pivoting( T,rr,cc )  ;
            end
        end
        Tu = T_new ;
        T = T_new ;
        
        %------------------------------------------------------------------
        %Da questo momento i valori di lambda e del
        %moltiplicatore di carico saranno incrementi rispetto allo step
        %precendete. I valori delle variabili Phi sono invece i valor
        %totali.
        %Proseguo facendo entrare in base il lambda relativo all'ultimo
        %modo attivato prima dello scarico
        cc = CCC(ii+1) ;  %colonna elemento Pivot (nota in partenza)
        
        Alpha = T(:,2)./T(:,cc) ;
        if last_mod < 0
            p = min(Alpha(Alpha>0)) ;
%             p = max(Alpha(Alpha<0)) ;
        else
            p = max(Alpha(Alpha<0)) ;
        end
        
        
        rr = find(Alpha==p) ; %riga elemento Pivot

        if size(rr,1)>1
            rr = rr(1);
        end
        
        % Controlli sulla riga pivot
        %Se un modo PhiB vuole uscire dalla base, devo prima verificare che
        %il relativo modo PhiA sia fuori dalla base. Se così non fosse,
        %passo a selezionare una diversa riga pivot.
        while T(rr,1)<= -10000 %un modo PhiB vuole uscire dalla base
            for ff = 1 : size(T,1)
                if T(rr,1)== -T(ff,1)
                    check = 1 ; %il relativo modo PhiA è ancora in base
                    Alpha(rr)=0 ; %passo a selezionare una nuova riga
                elseif T(rr,1)== -T(1,cc)
                    check = 1 ; %il relativo modo PhiA vuole entrare in base
                    mm(rr)=0 ; %passo a selezionare una nuova riga
                end
            end
            if last_mod < 0
                p = min(Alpha(Alpha>0)) ;
%                 p = max(Alpha(Alpha<0)) ;
            else
                p = max(Alpha(Alpha<0)) ;
            end
            if isempty(p) == 1
                display 'errore rapporto>0 - unbounded solution'
                break
            end
            rr = find(Alpha==p) ;  %riga elemento Pivot
            if size(rr,1)>1
                rr = rr(1);
            end
            
            if check == 0
                break
            else
                check = 0; %resetto il check
            end
        end
        %------------------------------------------------------------------
        rr 
        cc
        [T_new] = pivoting( T,rr,cc )  ;
        T_uu = T_new ;
        %-------------------- FINE FUNZIONE UNLOADING ---------------------
        
        %Se allo step prededente c'è stata una trasformazione vuota
        %Aggiorno il valore di Lambda che era stato azzerato.
        if Tcheck(ii) == 1
            T_new(rr,2) = T_new(rr,2) + Lambda_b ;
        end
            
        %Nel costruire il vettore dei modi attivi non considero
        %eventuali lambda rientrate in base in caso di scarichi olonomi
        if  T_new(1,cc) >= 10000|| T_new(1,cc) <= -10000
            act_mod = [act_mod; T_new(1,cc)] ;
        end
        
        %Aggiorno i valori dei modi attivi
        act_mod = [] ;
        for qq = 3 : size(T_new,2)
            if T_new(1,qq) >= 10000 || T_new(1,qq) <= -10000
                act_mod = [act_mod; T_new(1,qq)] ;
            end
        end
        
        last_mod = T_new(1,cc) ;
        uu = find(act_mod==T_new(1,cc)) ;
        prev_mod = act_mod ;
        prev_mod(uu) = [] ;
               
        Dalpha = T_new((find((T_new(:,1))==200000000)),2) ;
        alpha_s = alpha ;
        alpha = alpha_s + Dalpha ;
        Lambda_s = Lambda ;
        
%         %Se alpha è diventanto negativo, vuol dire che il passo è stato
%         %troppo lungo e devo scalare i valori per far arrivare alpha a zero
% 
%         if alpha<0
%             ss = alpha/Dalpha ;
%             alpha = 0 ;
%             Tdue = T_new(:,2) ;
%             T_new(:,2) = T_new(:,2)* ss;
%         end
        
    %----------------------------------------------------------------------
    % Caso di carico
    %----------------------------------------------------------------------
    else 
        rr
        cc
        [T_new] = pivoting( T,rr,cc ) ;
        if ii == 1
            Tp = T ;
        else
            if Tcheck(ii) ~= 2
                Tp = T ; %Se c'è stata una trasformazione vuota il Tp non viene aggiornato
            end
        end
        clear T
        
        %Se allo step prededente c'è stata una trasformazione vuota
        %Aggiorno il valore di Lambda che era stato azzerato.
        if Tcheck(ii) == 1
            T_new(rr,2) = T_new(rr,2) + Lambda_b ;
        end
        
        %Aggiorno i valori dei modi attivi
        act_mod = [] ;
        for qq = 3 : size(T_new,2)
            if T_new(1,qq) >= 10000 || T_new(1,qq) <= -10000
                act_mod = [act_mod; T_new(1,qq)] ;
            end
        end

%         %Nel caso di trasformazione vuota aggiungo ai modi attivi
%         %e a prev_mod il modo PhiA che ho eliminato dal tableu.
%         %Mi servirà nel caso di scarico per costruire il tableau dei modi
%         %attivi al passo precedente.
%         if Tcheck(ii) == 2
%             act_mod = [act_mod; T_new(1,cc)*10000] ; 
%         end

        %Aggiorno last_mod e prev_mod solo se non c'è stata una
        %trasformazione vuota. In quel caso rimangono invariati.
        
        if Tcheck(ii) ~= 2
            last_mod = T_new(1,cc) ;
            uu = find(act_mod==T_new(1,cc)) ;
            prev_mod = act_mod ;
            prev_mod(uu) = [] ;
        end


        if ii == 1
            alpha_old = alpha0 ;
        else
            alpha_old = alpha ;
        end
        
        Dalpha = T_new((find((T_new(:,1))==200000000)),2) ;
        
        %alpha_s e Lambda_s sono diverse da zero solo se c'è stato uno
        %scarico agli step precedenti
        alpha = alpha_s + Dalpha ;
        if Tcheck(ii) == 2 ;
            display 'Carico decrescente - Trasformazione vuota'
        elseif alpha - alpha_old <= 0
            checkload = 1 ; %ramo discendente del carico
            display 'carico decrescente'
        else
            display 'carico crescente'
        end
        
%         %Se alpha è diventanto negativo, vuol dire che il passo è stato
%         %troppo lungo e devo scalare i valori per far arrivare alpha a zero
% 
%         if alpha<0
%             ss = ALPHA(ii)/Dalpha ;
%             alpha = 0 ;
%             Tdue = T_new(:,2) ;
%             T_new(:,2) = T_new(:,2)* ss;
%         end
       
        
        %Check per vedere se alpha esce dalla base
        if isempty(T_new((find((T_new(:,1))==200000000)),2))== 1
            % se alpha esce dalla base
            alpha = 0 ;
            display 'error'
        end
    end
    Acheck = 0 ;  %Azzero il controllo
    
    
    %----------------------------------------------------------------------
    %%%%%%% CHECK SU TERMINI NOTI e ALPHA %%%%%%%%%
    %----------------------------------------------------------------------
    for jj = 2 : size(T_new,1)
        if T_new(jj,1) > 0 && T_new(jj,1) ~= 200000000
            if T_new(jj,2) < 0
                display 'ERRORE - VIOLAZIONE VINCOLI'
                vincoli = 1 ;
                break
            end
        end
    end
    
    if alpha < 0
        display 'ALPHA < 0'
        vincoli = 1
    end       
    
    if vincoli == 1
        break 
    end
    
    %----------------------------------------------------------------------
    %----------- Calcolo i vettori ALPHA, PHI e LAMBDA --------------------
    %----------------------------------------------------------------------
    % Costruzione del vettore ALPHA
    if ii == 1
        ALPHA = [alpha0 alpha] ;
    else
        ALPHA = [ALPHA alpha] ;
    end
    
    
    % Costruzione del vettore PHI
    Phi = zeros(size(P)) ;
    for pp = 1 : size(P,1)
        for tt = 1 : size(T_new,1)-1 %modi presenti nel tableau
            if P(pp,1)== T_new(tt+1,1)
                Phi(pp)=T_new(tt+1,2) ;
            end
        end
    end
    
    if ii == 1
        PHI = [Phi0 Phi] ;
    else
        PHI = [PHI Phi] ;
    end
    clear ff
    clear tt
    % Costruzione del vettore LAMBDA
    DLambda = zeros(size(L')) ;
    for ll = 1 : size(L,2)
        for tt = 1 : size(T_new,1)-1
            if L(1,ll)== T_new(tt+1,1)
                DLambda(ll)=T_new(tt+1,2) ;
            end
        end
    end
    Lambda = Lambda_s + DLambda ;
    
    if ii == 1
        LAMBDA = [Lambda0 Lambda] ;
        
        
    else
        LAMBDA = [LAMBDA Lambda] ;
    end
    clear ll
    clear tt
    %Definiamo un vettore lambda considerando per i modi con softening 
    %solo un modo per ogni nodo. Quando il modo B si attiva, la lambda 
    %totale coincide con LambdaB e LambdaA è stato eliminato dal Tableau.
    %I lambda relativi ai modi plastici rimangono invariati
    kk = (size(Lambda,1) - size(mp,1))/2 ;
    LambdaS = zeros(kk + size(mp,1),1) ;
    for jj = 1 : kk
        if Lambda(2*jj) > tol
            LambdaS(jj) = Lambda(2*jj-1) + Lambda(2*jj) ;
        else 
            LambdaS(jj) = Lambda(2*jj-1) ;
        end
    end
    for jj = 2*kk+1 : size(Lambda,1)
        LambdaS(jj-kk) = Lambda(jj) ;
    end
    
    LAMBDAS = [LAMBDAS LambdaS];
   
    clear ll
    
    
    %----------------------------------------------------------------------       
    %-------------------------- CALCOLO SPOSTAMENTI -----------------------
    %----------------------------------------------------------------------
    
   % Contributo elastico --------------------------------------------------
   
    %Calcolo gli spostamenti elastici
    Fc = alpha*F ;
    [u_c,R_c]=solve(K,Fc,bc);
    
    %Calcolo le forze nodali negli elementi di influenza a partire dagli
    %spostamenti elastici
    [Fnc] = F_nodal(Cpp,Bee,Be,C,et,exL,eyL,t,b,D,ndof,nne,ne,u_c) ;
       
    %Calcolo gli sforzi elastici
    es_c = []; %element stress matrix
    for i = 1 : ne
        ed_c(i,1:2*nne) = u_c(Dof(i,:))' ; %transform the nodal disp vector into element disp matrix
        [ess]=stress(et,ex(i,:),ey(i,:),D,ed_c(i,:)') ; %compute stress
        es_c = [es_c; ess] ;
    end
    clear ess i
    
    % Contributo dovuto alle lambda plastiche -----------------------------

    %Calcolo le deformazioni plastiche
    def = zeros(ne,3) ;
    for i = 1 : size(mp,1)
        y = mp(i) ;
        j = 2*size(Cp,1) + i ;
        deff = LambdaS(j) * NNN(i,:) ;
        def(y,:) = deff ;
    end
    clear i j y
    
    %Calcolo gli sforzi nella struttura a nodi fissi dovuti a lambdaS
    es_f_p = [] ;
    ef_f_p = [] ;
    for i = 1 : ne
        es_ff = -D*def(i,:)'  ; %Sforzi a nodi fissi
        es_f_p = [es_f_p; es_ff'] ; %Matrice degli sforzi della struttura
        [ef_ff] = fint(et,ex(i,:),ey(i,:),t,es_ff') ; %Forze nodali a nodi fissi
        ef_f_p = [ef_f_p; ef_ff] ; %Matrice delle forze nodali
    end
    clear i
    
    %Calcolo le forze nella struttura a nodi fissi dovute a lambdaS
    Fnf_p = ef_f_p' ;
    
    %Calcolo le forze nodali Fint (assemblate in un vettore globale),
    %cambio segno e applico alla struttura a nodi spostabili.
    Fint = zeros(ndof*nn,1) ;
    for i = 1 : ne
        for j = 1 : size(C,2)
            k = C(i,j) ;
            Fint(2*k-1) = Fint(2*k-1) + ef_f_p(i,2*j-1);
            Fint(2*k) = Fint(2*k) + ef_f_p(i,2*j);
        end
    end
    clear i j
    
    %Contributo spostamento e sforzo a nodi spostabili
    [u_s, R_s]=solve(K,-Fint,bc) ;
    es_s_p = []; %element stress matrix
    for i = 1 : ne
        ed_s(i,1:2*nne) = u_s(Dof(i,:))' ;
        [ess]=stress(et,ex(i,:),ey(i,:),D,ed_s(i,:)') ;
        es_s_p = [es_s_p; ess] ;
    end
    clear ess i
    
    % Calcolo le forze nodali a nodi spostabili elemento per elemento
    Fns_p = zeros(nne*ndof, ne) ; %ogni colonna mi da le forze nodali di ogni elemento
    for bbb = 1 : size(Cpp,1)
        for dd = Bee(bbb,1) : Bee(bbb,2)
            CC = C(Be(dd),:) ; %nodi dell'elemento considerato
            [Ke] = stif2d(et,exL(dd,:),eyL(dd,:),t,b,D) ;
            d = zeros(ndof*nne,1) ;
            for iii = 1 : nne %Number of element nodes
                for jj = 1 : ndof %Number of node dof
                    rw = ndof*(CC(iii)-1)+jj ; %Row position
                    d(ndof*(iii-1)+jj)= u_s(rw); %Displacement of element nodes
                end
            end
            clear iii jj bbb
            ff = Ke*d ;  %Forze sui nodi dell'elemento considerato
            ww = Be(dd) ;
            Fns_p(:,ww) = ff ;
        end
    end
    
    
    % Contributo dovuto alle lambda delle fessure -------------------------
      
    %Calcolo le forze Ft generate da Lambda sulla struttura a nodi fissi
    LambdaS_p = LambdaS(1:kk) ;
    [Ft, Fnf, df] = fixednodes(LambdaS_p,exL,eyL,C,Cp,Cpp,Be,Bee,N,NN,et,D) ;
    
    %Calcolo gli sforzi negli elementi dovuti alle distorsioni imposte
    es_f = []; %element stress matrix
    ed_f = df' ;
    for j = 1 : ne
        [ess]=stress(et,ex(j,:),ey(j,:),D,ed_f(j,:)') ; %compute stress
        es_f = [es_f; ess] ; %stress
    end
    clear j
    
    %Calcolo la soluzione elastica della struttura soggetta a -Ft
    [u_p,R_p]=solve(K,-Ft,bc) ;
    
    %Calcolo le forze nodali negli elementi di influenza a partire dagli
    %spostamenti
    [Fns] = F_nodal(Cpp,Bee,Be,C,et,exL,eyL,t,b,D,ndof,nne,ne,u_p) ;
    
%     %Proietto le forze tenendo conto degli elementi di influenza al nodo
%     [Fnp] = F_projection(F_np,ndof,nne,NN,Cp,Cpp,Be,Bee,C) ;
    
    %Calcolo gli sforzi negli elementi dovuti alle distorsioni imposte
    es_s = []; %element stress matrix
    for j = 1 : ne
        ed_s(j,1:2*nne) = u_p(Dof(j,:))' ; %transform the nodal disp vector into element disp matrix
    end
    clear j
    for j = 1 : ne
        [ess]=stress(et,ex(j,:),ey(j,:),D,ed_s(j,:)') ; %compute stress
        es_s = [es_s; ess] ; %stress
    end
    clear j
    
    %Proiezione delle lambda su x e y
    %!!!!!!!!!1 I valori di spostamenti totale per alcuni elementi sono
    %SBAGLIATI (per gli elementi comuni a due linee di fessura diverse)!!!
    
    
    for ccc = 1 : size(Cp,1)
        RRR = NN(2*ccc-1:2*ccc,:) ;
        RRRR = RRR(:,1)' ;
        Lambdax(ccc,1) = RRRR * LambdaS(2*ccc-1:2*ccc) ;
        RRRR = RRR(:,2)' ;
        Lambday(ccc,1) = RRRR * LambdaS(2*ccc-1:2*ccc) ;
    end
    %Spostamenti non compatibili
    u_t = zeros(nn*ndof,1) ;
    u_t(2*Cp-1) = Lambdax ;
    u_t(2*Cp) = Lambday ;
    
    
       
%         %Verifica
%         [bF] = vector_b(et,ex,ey,N,C,Fc,b,bc,D,Cp,Cpp,Be,Bee,NN) ;
%         bbF = zeros(2*size(bF,1),1) ;
%         for xx = 1 : size(bF,1)
%             bbF(2*xx-1) = bF(xx);
%             bbF(2*xx) = bF(xx); 
%         end
%         clear xx
%         
%         V = -A*Lambda+bbF ;        
%         if ii == 1
%             VV = [V] ;
%         else
%             VV = [VV V] ;
%         end


    % Spostamenti totali --------------------------------------------------  
   
    u = u_c + u_s + u_p ; %compatibili
    
    U_TOT = [U_TOT u] ;
    
    %Per plottare la deformata bisogna assegnare lo spostamento dovuto al
    %Lambda della fessurazione solo ad un lato della struttura
    
    for i = 1 : ne
        ed(i,1:2*nne) = u(Dof(i,:))' ; %transform the nodal disp vector into element disp matrix
    end
    
    d = zeros(size(ed));
    for ccc = 1 : size(Cp,1)
        for bbb = 1 : size(Cpp,1)
            if (ccc >= Cpp(bbb,1)) && (ccc <= Cpp(bbb,2))
                for dd = Bee(bbb,1) : Bee(bbb,2)
                    for iii = 1 : nne
                        CC = C(Be(dd),iii) ;   % nodo elemento considerato
                        ww = Be(dd) ;
                        if Cp(ccc) == CC
                            d(ww,2*iii-1) = d(ww,2*iii-1)  - 20* Lambdax(ccc) ; %amplifico lambda per il plot
                            d(ww,2*iii) = d(ww,2*iii) - 20*Lambday(ccc)  ;
                        end
                    end
                end
            end
        end
    end
    ed = ed + d ;
    clear ccc bbb dd iii ww
    
    
%     if ii==1 || ii> 60
%         k=2;
%         h = figure ;
%         x=ex' ;
%         y=ey';
%         xc=[x; x(1,:)];
%         yc=[y; y(1,:)];
%         x_d=(ex+k*ed(:,[1 3 5]))';
%         y_d=(ey+k*ed(:,[2 4 6]))';
%         xc_d=[x_d; x_d(1,:)];
%         yc_d=[y_d; y_d(1,:)];
%         hold on
%         axis equal
%         plot(xc,yc,'Color',[0.7 0.7 0.7]);
%         plot(xc_d,yc_d,'k') ;
%         set(0,'defaultlinelinewidth',1.5);
%         ax.GridLineStyle = '-';
%         % grid on
%         box on
%         % xlim([-200 1400]) ;
%         % ylim([-200 1400]) ;
%         saveas(h,sprintf('FIG_DEf%d.png',ii));
%     end
    
    
    % Forze totali ---------------------------------------------------------
    
    F_tot = Fnc + Fnf_p + Fns_p + Fnf + Fns ;
    
    %Proietto le forze tenendo conto degli elementi di influenza al nodo
    [Fn] = F_projection(F_tot,ndof,nne,NN,Cp,Cpp,Be,Bee,C) ;
        
    FN = [FN Fn] ;
    
    Rf = R_c + R_s + R_p ;
    
    RF = [RF Rf] ;
    
    es_tot = es_c + es_f_p + es_s_p + es_f + es_s ;
    
    %Sforzi principali
    es_1 = ((es_tot(:,1)+es_tot(:,2))/2) + sqrt(((es_tot(:,1)-es_tot(:,2))/2).^2 + es_tot(:,3).^2) ;
    es_2 = ((es_tot(:,1)+es_tot(:,2))/2) - sqrt(((es_tot(:,1)-es_tot(:,2))/2).^2 + es_tot(:,3).^2) ;
    
    
    %Proietto gli sforzi
    es_totP = [] ;
    for j = 1 : size(mp,1)
        i = mp(j) ;
        es_totp = es_tot(i,:) * NNN(j,:)' ;
        es_totP = [es_totP; es_totp];
    end
    clear j i
    
    ES_TOT = [ES_TOT es_tot] ;
    ES_TOTP =[ES_TOTP es_totP] ;
 
% % %    Per il plot del meccanismo di collasso     
% %     if T(:,cc)>=-tol %|T(find(T(:,1)==200000000),cc)<tol
% %         display 'collasso'
% %         Lambda_mec = zeros(size(Lambda)) ;
% %         Lambda_mec(T(1,cc)) = 1 ;
% %         for ll = 1 : size(L,2)
% %             for tt = 1 : size(L,2)
% %                 if L(1,ll)== T(tt+1,1)
% %                     Lambda_mec(ll)=T(tt+1,cc) ;
% %                 end
% %             end
% %         end
% %         [ed,u,u_t] = displacement(alpha,F,K,bc,Lambda_mec,exL,eyL,C,Cp,Cpp,Be,Bee,N,NN,et,D,Dof) ;    
% %         break
% %     end
%              
 
    if alpha==0
        break
    end


end



%--------------------------------------------------------------------------        
%-------------------------- PLOT DIAGRAMMI --------------------------------
%-------------------------------------------------------------------------- 
% Costruisco il vettore ETT per eliminare dai grafici gli step delle 
% trasformazioni vuote in cui i valori di FN non sono corretti. Nel caso del
% vettore dello spostamento e di Alpha uso ETT+1 perchè il vettore ha una 
% dimesione maggiore visto che considero anche lo step elastico.
ETT = 0 ;
for xx = 1 : size(ET,1)
    if ET(xx) == 1
        if ETT == 0
            ETT = xx ;
        else
            ETT = [ETT xx] ;
        end
    end
end
clear xx

% % Apertura di fessura vs Forza normale ------------------------------------
% XX = LAMBDAS(1:2*size(Cp,1),:) ;
% YY = FN ;
% 
% %Elimino dal plot gli step relativi alle trasformazioni vuote che non hanno
% %significato meccanico
% if ETT ~= 0 
%     XX(:,ETT) = XX(:,ETT-1);
%     YY(:,ETT+1) = YY(:,ETT) ;
% end
% 
% X = [zeros(size(XX,1),2) XX] ;
% Y = [zeros(size(YY,1),1) YY]/1000 ;
% 
% for i = 1 : size(X,1)
%     h = figure ;
%     hold on
%     plot(X(i,:),Y(i,:),'-o',...
%         'LineWidth',1.5,...
%         'Color','k',...
%         'MarkerEdgeColor','k',...
%         'MarkerFaceColor','b',...
%         'MarkerSize',5);
% %     for j = 1 : step
% %         text(X(i,j+2),Y(i,j+2),['   ' num2str(j)],'FontSize',10) ;
% %     end
%     ax.GridLineStyle = '-';
%     set(gca,'FontSize',17)
%     set(gca, 'FontName', 'Times New Roman')
%     %     title(['Node ', num2str(Cp(i))]);
%     if rem(i,2)==0
%         xlabel('Crack sliding \lambdas [mm]', 'FontSize',22);
%         ylabel('Force Ft [kN]','FontSize',22);
%     else
%         xlabel('Crack opening \lambdaw [mm]', 'FontSize',22);
%         ylabel('Force Fn [kN]','FontSize',22);
%     end
%     
%     %     xlim([0 1]) ;
%     %     ylim([0 30000]) ;
%     grid on
%     box on
%     saveas(h,sprintf('FIG_FLambda%d.png',i));
%     hold on
% end
% 
% 
% % Lambda vs Sforzo ------------------------------------
% X = LAMBDAS(2*size(Cp,1)+1:size(LAMBDAS,1),:) ;
% Y = ES_TOTP ;
% 
% %Elimino dal plot gli step relativi alle trasformazioni vuote che non hanno
% %significato meccanico
% if ETT ~= 0 
%     X(:,ETT) = X(:,ETT-1) ;
%     Y(:,ETT+1) = Y(:,ETT) ;
% end
% 
% X = [zeros(size(X,1),2) X] ;
% Y = [zeros(size(Y,1),1) Y] ;
% for i = 1 : size(X,1)
%     hh = figure ;
%     hold on
%     plot(X(i,:),Y(i,:),'-o',...
%                 'LineWidth',1.5,...
%                 'Color','k',...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','b',...
%                 'MarkerSize',5);
%     for j = 1 : step
%         text(X(i,j+2),Y(i,j+2),['   ' num2str(j)],'FontSize',10) ;
%     end
%     ax.GridLineStyle = '-';
%         set(gca,'FontSize',17)
%     set(gca, 'FontName', 'Times New Roman')
% %     title(['Element ', num2str()]);
%     set(gca, 'FontName', 'Times New Roman')
%     xlabel('Plastic deformation \lambdap [mm]', 'FontSize',22);
%     ylabel('Resistance \sigma*n [MPa]','FontSize',22);
% %     xlim([0 1]) ;
% 
% 
% %     ylim([0 30000]) ;
%     grid on
%     box on
%     saveas(hh,sprintf('FIG_sigmalambda%d.png',i));
%     hold on
% end

% 
% Carico vs Spostamento --------------------------------------------
load('Experimental.txt');
Exp = Experimental ;
x_e = Exp(:,1) ;
y_e = Exp(:,2);

dd = 10000 ; %numero figura
% xx = abs(U_TOT(26,:)) + abs(U_TOT(58,:)) ;
xx = abs( U_TOT(34,:) - U_TOT(18,:))   ;
% xx = abs( U_TOT(58,:) - U_TOT(26,:))   ;
% xx = abs( U_TOT(138,:) - U_TOT(26,:))   ;
% xx = (abs(U_TOT(58,:) - U_TOT(26,:)) + abs(U_TOT(62,:) - U_TOT(22,:)))/2  ;
% xx = (abs(U_TOT(138,:) - U_TOT(26,:)) + abs(U_TOT(78,:) - U_TOT(86,:)))/2  ;
% xx = (U_TOT(14,:) + U_TOT(16,:) + U_TOT(18,:))/3   ;
% xx = (U_TOT(42,:) + U_TOT(44,:) + U_TOT(46,:)+ U_TOT(48,:)+ U_TOT(50,:))/5   ;
yy = ALPHA ;
if ETT ~= 0
    xx(ETT+1) = xx(ETT);
    yy(ETT+1) = yy(ETT) ;
end
xx = [0 xx]' ;
yy = [0 yy]' ;
figure(dd)
plot(xx,yy,'-o',...
                'LineWidth',1.5,...
                'Color','k',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',4);
% for j = 1 : step  
%     text(xx(j+2),yy(j+2),['   ' num2str(j)],'FontSize',8) ;
% end
ax.GridLineStyle = '-';
% title('Node 25');
xlabel('Dispacement [mm]','FontSize',14);
ylabel('Load [N]','FontSize',14);
grid on
box on
hold on
plot(x_e,y_e,'--') ;


% saveas(gcf,'FIG_loadvsdisp.png')

% % % Moltiplicatore vs Step --------------------------------------------
% dd = 10000 ; %numero figura
% 
% xx = [1 : step + 1] ;
% yy = ALPHA ;
% if ETT ~= 0
% %     xx(ETT+1) = xx(ETT);
%     yy(ETT+1) = yy(ETT) ;
% end
% xx = [0 xx]' ;
% yy = [0 yy]' ;
% figure(dd)
% plot(xx,yy,'-o',...
%                 'LineWidth',1.5,...
%                 'Color','k',...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','b',...
%                 'MarkerSize',4);
% % for j = 1 : step  
% %     text(xx(j+2),yy(j+2),['   ' num2str(j)],'FontSize',7) ;
% % end
% ax.GridLineStyle = '-';
% % title('Node 25');
% xlabel('Dispacement [mm]','FontSize',14);
% ylabel('Load [N]','FontSize',14);
% grid on
% box on
% saveas(gcf,'FIG_loadvsdisp.png')

% Deformata --------------------------------------------------------------
figure('Name','Deformed shape')
k=2;
x=ex' ;
y=ey';
xc=[x; x(1,:)];
yc=[y; y(1,:)];
x_d=(ex+k*ed(:,[1 3 5]))';   
y_d=(ey+k*ed(:,[2 4 6]))';
xc_d=[x_d; x_d(1,:)];
yc_d=[y_d; y_d(1,:)];
hold on
axis equal
plot(xc,yc,'Color',[0.7 0.7 0.7]);
plot(xc_d,yc_d,'k') ;
set(0,'defaultlinelinewidth',1.5);
ax.GridLineStyle = '-';
% grid on
box on
% xlim([-200 1400]) ;
% ylim([-200 1400]) ;
saveas(gcf,'figure18.png')

% yticks(-1000:200:1000) ;

% colormap('jet')
% fill(ex',ey',es_1')
% axis equal
toc

