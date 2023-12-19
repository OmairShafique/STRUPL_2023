function [T_A] = Act_mod( prev_mod, last_mod, A, L_p, bb )
%La funzione assembla il tableau dei modi attivi

prev_mod = sort(prev_mod) ;

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

% Assemblo il tableau
T_A = [0 100000000 prev_mod'/10000 last_mod/10000; prev_mod -bp Ap apl; last_mod -bl apl' all] ; 


end

