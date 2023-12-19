function [T_La,T_L] = ACheck( act_mod,DAlpha,T,Tp)
%ACHECK costruisce il tableau dei modi attivi e controlla che la matrice A
%sia semi-definita positiva

last_mod = act_mod(end) ;  %ultimo modo attivato
prev_mod = act_mod(1:(end-1)) ; %altri modi precedentemente attivati
prev_mod = sort(prev_mod) ;

%Sostituisco al valore di alpha dell'ultimo tableau il DAlpha o zero ed effettuo una
%trasformazione facendo entrare in base l'ultimo modo attivato e
%facendo uscire dalla base il DAlpha (=scrivo il tableau all'inizio
%dello step ii)
rr_n = (find((T(:,1))==20000)) ;
cc_n = (find((T(1,:))==last_mod)) ;
T_newL = T ;
T_newL(rr_n,2) = 0  ;
[T_L] = pivoting(T_newL,rr_n,cc_n)  ; % tableau con Alpha fuori base

%Calcolo il vettore dei termini noti (incremento) del tableau alla Lemke come t(I)-t(I-1)/DAlpha
%Termini noti nel tableau allo step I-1
tn_p = [zeros(size(prev_mod)); 0] ;
for jj = 1 : size(Tp,1)
    for kk = 1 : size(prev_mod,1)
        if Tp(jj,1) == prev_mod(kk,1)/100
            tn_p(kk) = Tp(jj,2) ;
        end
    end
end
clear kk jj
for  jj = 1 : size(Tp,1)
    if Tp(jj,1) == last_mod
        tn_p(end) = Tp(jj,2) ;
    end
end
%Termini noti nel tableau allo step I
tn_l = [zeros(size(prev_mod)); 0] ;
for jj = 1 : size(T,1)
    for kk = 1 : size(prev_mod,1)
        if T(jj,1) == prev_mod(kk,1)/100
            tn_l(kk) = T(jj,2) ;
        end
    end
end
clear kk jj
for  jj = 1 : size(T,1)
    if T(jj,1) == last_mod
        tn_l(end) = T(jj,2) ;
    end
end

tn = (tn_l-tn_p)/DAlpha;

% Trovo la posizione dei coefficenti relativi ai modi attivati nella matrice Tp
% e individuo i valori da posizionare nel tableau alla Lemke
prev_mod(find(prev_mod==0)) = [] ;
for pp = 1 : size(prev_mod,1)
    for aa = 1 : size(T_L,2)
        if prev_mod(pp) == T_L(1,aa)
            p_c(pp,1) = aa ;
        end
    end
end
clear pp aa

for pp = 1 : size(prev_mod,1)
    for aa = 1 : size(T_L,1)
        if prev_mod(pp)/100 == T_L(aa,1)
            p_r(pp,1) = aa ;
        end
    end
end
clear pp aa

for pp = 1 : size(last_mod,1)
    for aa = 1 : size(T_L,2)
        if last_mod(pp)/100 == T_L(1,aa)
            l_c(pp,1) = aa ;
        end
    end
end
clear pp aa

for pp = 1 : size(last_mod,1)
    for aa = 1 : size(T_L,1)
        if last_mod(pp) == T_L(aa,1)
            l_r(pp,1) = aa ;
        end
    end
end
clear pp aa

L1 = T_L(p_r, p_c) ;
L2 = T_L(l_r, l_c) ;
L3 = T_L(p_r, l_c) ;
L4 = T_L(l_r, p_c) ;

% Assemblo il tableau
T_La = [0 10000 prev_mod' last_mod/100; prev_mod/100 tn(1:(end-1)) L1 L3; last_mod tn(end) L4 L2] ;

end

