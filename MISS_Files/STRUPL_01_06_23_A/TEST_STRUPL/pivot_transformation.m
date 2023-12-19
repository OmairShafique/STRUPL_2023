function [T_new] = pivot_transformation( T,rr,cc )
%PIVOTING La funzione effettua una trasformazione pivotale
%   La funzione, dato un tableau T di partenza e date la posizione della
%   riga pivot rr e della colonna pivot tt, restituisce il nuovo tableau, 
%   calcolandone i singoli elementi secondo le note regole.

T_new = zeros(size(T));
T_new(1,:) = T(1,:) ; 
T_new(:,1) = T(:,1) ; 
T_new(1,cc)= T(rr,1) ; %Elemento che entra in base
T_new(rr,1)= T(1,cc) ; %Elemento che esce dalla base

 for i = 2 : size(T,1)  %elemento generico
     for j = 2 : size(T,2) 
        T_new(i,j) = T(i,j) - 1/T(rr,cc) * T(i,cc) * T(rr,j) ;
     end
 end
 
 clear i 
 clear j
 
 for j = 2 : size(T,2)  %riga elemento Pivot
     T_new(rr,j) = - 1/T(rr,cc) * T(rr,j) ;    
 end
 
 clear j

 for i = 2 : size(T,1)  %colonna elemento Pivot
     T_new(i,cc) = 1/T(rr,cc) * T(i,cc)  ; 
 end
 
 clear i
 
 T_new(rr,cc) = 1/T(rr,cc)  ;

end

load 'ragnettodata.mat' ,A,BI,R
% initialization

%Rappresentazione della tabella iniziale

A1= A(:,1);
A2= A(:,2);
A3= A(:,3);
A4= A(:,4);
A5= A(:,5);
T_1= table(R,-BI,A1,A2,A3,A4,A5,'VariableNames',{'R','α','λ1+','λ1-','λ3+','λ3-','λ6+'}, 'RowNames',{'φ1+','φ1-','φ3+','φ3-','φ6+'})
T=table2array(T_1);

% Ricerca dell'elemento Pivot
Alfa =- T(:,1)./T(:,2) ;
p = min(Alfa(Alfa>0)) ;
cc = 2 ; %colonna elemento Pivot (nota in partenza)
rr = find(Alfa==p) ; %riga elemento Pivot
if size(rr,1)>1
    rr = rr(1);
end

% Costruzione del nuovo tableau T_new -------------------------------------

% Il nuovo tableau si costruisce applicando le regole di pivoting
% Il primo elemento ad entrare in base è sempre Alfa
% Si utilizza la funzione 'pivoting' che effettua la trasformazione
% restituendo il nuovo tableau

[T_new] = pivoting( T,rr,cc )
T_new_t = array2table(T_new,'VariableNames',{' ','φ1+','λ1+','λ1-','λ3+','λ3-','λ6+'}, 'RowNames',{'α','φ1-','φ3+','φ3-','φ6+'})
T_new=table2array(T_new_t);
Alfa0 = T_new((find((T_new(:,1))==Alfa)),1) ;
% Ricerca dell'elemento Pivot
Alfa2 =- T_new(:,1)./T_new(:,3) ;
p1 = min(Alfa2(Alfa2>0)) ;
cc1 = 3 ; %colonna elemento Pivot (nota in partenza)
rr1 = find(Alfa2==p1) ; %riga elemento Pivot
if size(rr1,1)>1
    rr1 = rr1(1);
end

% Costruzione del nuovo tableau T_2 -------------------------------------

% Il nuovo tableau si costruisce applicando le regole di pivoting
% Il primo elemento ad entrare in base è sempre Alpha
% Si utilizza la funzione 'pivoting' che effettua la trasformazione
% restituendo il nuovo tableau

[T_2] = pivoting( T_new,rr1,cc1 )
T_2_t = array2table(T_2,'VariableNames',{' ','φ1+','φ3+','λ1-','λ3+','λ3-','λ6+'}, 'RowNames',{'α','φ1-','λ1+','φ3-','φ6+'})
T_2=table2array(T_2_t);
Alfa1 = T_2((find((T_2(:,1))==Alfa2)),3) ;
