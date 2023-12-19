function [T_A_nh] = Act_mod_no_olo(T_A)
%La funzione assembla il tableau dei modi non olonomi attivi

nh = [] ;
h = [] ;
for i = 2 : size(T_A,1)
    if T_A(i,1) > 0
        nh = [nh; i] ;
    elseif T_A(i,1) < 0
        h = [h; i] ;
    end
end


A_nh_nh = T_A(nh,nh+1) ;
A_h_h = T_A(h,h+1) ;
A_nh_h = T_A(nh,h+1) ;
A_h_nh = T_A(h,nh+1) ;
b_nh = T_A(nh,2) ;
b_h = T_A(h,2) ;

A = A_nh_nh - A_nh_h * A_h_h^-1 * A_h_nh ;
b = A_nh_h * A_h_h^-1 * b_h - b_nh ;
 
% Assemblo il tableau
T_A_nh = [0 100000000 T_A(nh,1)'/10000; T_A(nh,1) b A] ; 


end

