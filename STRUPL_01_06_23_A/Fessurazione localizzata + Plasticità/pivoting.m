function [T_new] = pivoting( T,rr,cc )
%PIVOTING La funzione effettua una trasformazione pivotale
%   La funzione, dato un tableau T di partenza e date la posizione della
%   riga pivot rr e della colonna pivot tt, restituisce il nuovo tableau, 
%   calcolandone i singoli elementi secondo le note regole.

T_new = zeros(size(T));
T_new(1,:) = T(1,:) ; 
T_new(:,1) = T(:,1) ; 
T_new(1,cc)= T(rr,1) ; %Elemento che entra in base
T_new(rr,1)= T(1,cc) ; %Elemento che esce dalla base

 for i = 2 : size(T,1) , %elemento generico
     for j = 2 : size(T,2) , 
        T_new(i,j) = T(i,j) - 1/T(rr,cc) * T(i,cc) * T(rr,j) ;
     end
 end
 
 clear i 
 clear j
 
 for j = 2 : size(T,2) , %riga elemento Pivot
     T_new(rr,j) = - 1/T(rr,cc) * T(rr,j) ;    
 end
 
 clear j

 for i = 2 : size(T,1) , %colonna elemento Pivot
     T_new(i,cc) = 1/T(rr,cc) * T(i,cc)  ; 
 end
 
 clear i
 
 T_new(rr,cc) = 1/T(rr,cc)  ;

end

