function [pos, AAA, v,AAA_v,AAA_b ] = def_pos_test_fes( prev_mod, last_mod, A,bb )
%La funzione controlla la positività della matrice dei modi attivi
%pos = 2 indica un derminante negativo
%pos = 1 indica un derminante = 0
%pos = 0 indica un determinante >0

prev_mod = sort(prev_mod) ;

for i = 1 : size(prev_mod,1)
    if  prev_mod(i)/10000 > 0
        pm(i,1) = 2*prev_mod(i)/10000 - 1;
    elseif prev_mod(i)/10000 < 0
        pm(i,1) = -(2*prev_mod(i)/10000);
    end
end

if last_mod/10000 < 0
    lm = -2*last_mod/10000;
elseif  last_mod/10000 > 0
    lm = 2*last_mod/10000 - 1;
end

Ap = A(pm,pm) ;
apl = A(pm,lm) ;
all = A(lm,lm) ;
bp = bb(pm) ;
bl = bb(lm) ;

AAA = [Ap apl; apl' all];
v = [-Ap^(-1)*apl;1] ;
AAA_v = AAA * v ;
AAA_b = -apl'*Ap^(-1)*bp + bl ;

pos = 0 ;

for i = 1 : size(AAA,1)
    if det(AAA(1:i, 1:i)) == 0
        pos = 1 ;
    elseif det(AAA(1:i, 1:i)) < 0
        pos = 2 ;
    end
end

end

