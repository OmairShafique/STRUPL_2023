function [T_La,Ap] = ACheck_pl_fes( prev_mod, last_mod, L_p, A,bb)
%ACHECK costruisce il tableau dei modi attivi e controlla che la matrice A
%sia semi-definita positiva
% %Vettore dei modi attivi diviso in ultimo modo attivato e modi precedentemente attivati
% last_mod = act_mod(size(act_mod,1),1)   %ultimo modo attivato
% prev_mod = act_mod(1:(size(act_mod,1)-1),1) ; %altri modi precedentemente attivati


prev_mod = sort(prev_mod) ;

% prev_mod
last_mod

for i = 1 : size(prev_mod,1)
    if prev_mod(i)/10000 < L_p(1) && prev_mod(i)/10000 > 0
        pm(i,1) = 2*prev_mod(i)/10000 - 1;
    elseif prev_mod(i)/10000 < 0
        pm(i,1) = -(2*prev_mod(i)/10000);
    elseif prev_mod(i)/10000 >= L_p(1)
        pm(i,1) = prev_mod(i)/10000 + L_p(1) -1 ;
    end
end

if last_mod/10000 < L_p(1) && last_mod/10000 < 0
    lm = -2*last_mod/10000;
elseif last_mod/10000 < L_p(1) && last_mod/10000 > 0
    lm = 2*last_mod/10000 - 1;
elseif last_mod/10000 >= L_p(1)
    lm = last_mod/10000 + L_p(1) -1 ;
end

Ap = A(pm,pm) ;
apl = A(pm,lm) ;
all = A(lm,lm) ;
bp = bb(pm) ;
bl = bb(lm) ;

det(Ap)
tnp = Ap^(-1)*bp ;
tnl = -(bl-apl'*Ap^(-1)*bp) ;
pp = Ap^(-1) ;
lp = apl'*Ap^(-1) ;
pl =  -Ap^(-1)*apl ;
ll = all - apl'*Ap^(-1)*apl;

% Assemblo il tableau
T_La = [0 100000000 prev_mod' last_mod/10000; prev_mod/10000 tnp pp pl; last_mod tnl lp ll] ; 

end

