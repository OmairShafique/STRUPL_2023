%%% Sdoppiare le righe per considerare due elmenti per ogni nodo

clc 
clear all

load('Influence_lenght_1.txt'); %lunghezze di influenza ai nodi
ILL = Influence_lenght_1;
IL = zeros(2*size(ILL,1),2) ;
IL(1:2:end)= ILL(1:end) ;
IL(2:2:end)= ILL(1:end) ;