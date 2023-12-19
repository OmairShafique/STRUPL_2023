function [T_new] = pivoting( T,rr,cc )
%PIVOTING La funzione effettua una trasformazione pivotale
%   La funzione, dato un tableau T di partenza e date la posizione della
%   riga pivot rr e della colonna pivot tt, restituisce il nuovo tableau, 
%   calcolandone i singoli elementi secondo le note regole.

T_new = zeros(size(T));
T_new(1,:) = T(1,:); 
T_new(:,1) = T(:,1);  
T_new(1,cc)= T(rr,1);  %Elemento che entra in base
T_new(rr,1)= T(1,cc);  %Elemento che esce dalla base
t_0_i = T(rr,:);
t_0_j = T(:,cc)
t_0 = t_0_j .* t_0_i

 for i = 1 : size(T,1)  %elementi generici
     for j = 1 : size(T,2)  
        T_new = T - 1/T(rr,cc) .* t_0 ;
     end
 end
 
 clear i 
 clear j
 
 for i = 1 : size(T,1) , %riga elemento Pivot
     t_0_i = T(rr,:);
     t_new_i = - 1/T(rr,cc) .* t_0_i ;    
 end
 
 clear i

 for j = 1 : size(T,2) , %colonna elemento Pivot
     t_0_j = T(:,cc);
     t_new_j = 1/T(rr,cc) .* t_0_j ; 
 end
 
 clear i
T_new(rr,:)= t_new_i
T_new(:,cc)= t_new_j
T_new(rr,cc) = 1/T(rr,cc)  

end




