% Programma per il calcolo del Ragnetto considerando indipendentemente ogni 
% variabile olonoma o non olonoma mediante uno schema tipo PLCP olonomo; 
% le variabili sono le Lambda e le Fi totali, variabili tutte positive.
% La soluzione iniziale prevede tutte le variabili Fi in base e tutte le
% variabili Lamda fuori base insieme al parametro ALFA.
% Si suppone che tutte le Resistenze R siano tutte positive (è facile
% estendere poi al caso in cui qualche Resistenza sia uguale a zero).
% Nel primo step si porta in base il parametro ALFA (colonna j) e si fa uscire 
% dalla base la Fi che per prima si annulla (riga i).
% Nello Step generico una variabile Lambda (i) entra in base e la
% variabile coniugata Fi(j) esce dalla base.
% Quando una o più variabili non olonome vogliono scaricarsi,cioè quando 
% qualche termine della colonna pivot corrispondente alle variabili Lambda
% non olonome sono negativi,si risolve un problema di LCP del problema 
% incrementale. L'idea è quella di risolverlo partendo dal Tableau
% iniziale, cioè quello in cui la colonna corrispondente alla variabile
% fuori base che dovrebbe entrare in base con alcuni(almeno uno)termini 
% corrispondenti a variabili Lambda non olonome in base negativi. 
% A tale Tableau si aggiunge una colonna tutta positiva corrispondente al nuovo 
% parametro E, seguento uno schema tipo Lemke. 
% Sarebbe opportuno a questo punto ricordare le ragioni meccaniche per 
% l'utilizzo di tale schema.
% Il procedimento dovrebbe prevedere:
% (i) si fa entrare in base il parametro E e si fa uscire il parametro ALFA.
% A tale scopo si dovrà prevedere la colonna aggiuntiva tutta positiva con
% un termine negativo opportuno corrispondente alla variabile ALFA.
% Dopo aver fatto entrare E e uscire ALFA si procede facendo entrare la
% variabile Lambda(j) che per ultima si era attivata nel processo di carico
% e aveva fatto si che fuori base si trovassero contemoraneamente le
% variabili Fi e Lambda coniugate. 
% Tale Lambda deve rimanere in base.
% Il processo continua come sempre controllando solo le Fi e le Lambda
% corrispondenti ai modi non olonomi attivi, che cioè hanno una Fi=0
% all'inizio della procedura di scarico.
% Il processo termina quando il parametro E esce dalla base.
% A questo punto si fa rientrare in base il parametro ALFA e il processo
% continua come previsto nella sua formulazione olonoma.
%
% Il presente programma viene particolarizzato al problema del ragnetto per
% controllare il funzionamento della procedura sia nel caso di scarico per
% carico crescente proporzionalmente, sia nel caso di formazione di
% meccanismo (collasso) e pseudo meccanismo (scarico).
%
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



