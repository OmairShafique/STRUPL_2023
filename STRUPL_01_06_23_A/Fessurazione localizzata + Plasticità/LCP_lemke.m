function [Tz_fin,Err, Tz_1] = LCP_lemke(T)
%Risolve il problema di LCP col metodo Lemke

Err = 0; 

%Aggiungo la colonna di z0 (tutti 1)
zo = ones(size(T,1)-1,1) ;
Tz = [T [50000000; zo]] ;
%Eseguo la prima trasformazione facendo entrare in base z0
cc = size(Tz,2) ;
mm= Tz(:,2)./Tz(:,cc) ;
mm(1)= 100000000 ;  %riga di intestazione
mm(end) = [] ; %non considero l'ultima riga

p = min(mm(mm<0)) ;
rr = find(mm==p) ;
if size(rr,1)>1
    rr = rr(1);
end

%Se un modo LambdaB vuole uscire dalla base,seleziono una diversa riga pivot.
check = 0 ;  %il check mi serve per interrompere il ciclo quando necessario
while Tz(rr,1)< 0 && Tz(rr,1)>-10000  %un modo LambdaB vuole uscire dalla base
    check = 1 ;
    display 'LCP - LambdaB esce dalla base'
    mm(rr)=10000000000 ; %passo a selezionare una nuova riga
    p = min(mm(mm<0)) ;
    if isempty(p) == 1
        display 'errore rapporto>0 - unbounded solution'
        break
    end
    rr = find(mm==p) ;  %riga elemento Pivot
    if size(rr,1)>1
        rr = rr(1);
    end
    if check == 0
        break
    else
        check = 0; %resetto il check
    end
end
    
[Tz_new] = pivoting( Tz,rr,cc ) ;


%Eseguo le trasformazioni successive fino a che z0 è fuori base
for gg = 1 : 10000

    Tz = Tz_new ;
    if Tz(1,cc) == 50000000
        Tz_fin = Tz ;
        Tz_fin(:,cc) = [] ;
        break
    elseif Tz(1,cc)>=10000 || Tz(1,cc)<=-10000
        cc = find(Tz(1,:)==Tz(1,cc)/10000) ;
    else
        cc = find(Tz(1,:)==Tz(1,cc)*10000) ;
    end
    mm= Tz(:,2)./Tz(:,cc) ;
    mm(1)= 10000000 ;  %riga di intestazione
    mm(end) = [] ; %ultima equazione di controllo non considerata
    p = max(mm(mm<0)) ;
    if isempty(p) == 1
        display 'LCP unbounded solution'
        Err = 1 ;
        Tz_fin = Tz ;
        break
    end
    rr = find(mm==p) ; %riga elemento Pivot
    if size(rr,1)>1
        rr = rr(1);
    end
    %Controllo eventuali zero nella colonna dei termini noti
    %Senza il check il rapporto tra zero ed un numero negativo
    %non viene considerato.
    for j =  1 : size(mm, 1)
        if Tz(j,2) == 0 && Tz(j,cc)<0
            rr = j ;
            break
        end
    end
    
    if Tz(rr,1)< 0 && Tz(rr,1)>-10000 
        display 'LCP - LambdaB esce dalla base'
    end
    
    %Se un modo LambdaB vuole uscire dalla base,seleziono una diversa riga pivot.
    check = 0 ;  %il check mi serve per interrompere il ciclo quando necessario
    while Tz(rr,1)< 0 && Tz(rr,1)>-10000  %un modo PhiB vuole uscire dalla base
        check = 1 ;
        mm(rr)=0 ; %passo a selezionare una nuova riga
        p = max(mm(mm<0)) ;
        if isempty(p) == 1
            display 'errore rapporto>0 - unbounded solution'
            Err = 1 ;
            break
        end
        rr = find(mm==p) ;  %riga elemento Pivot
        if size(rr,1)>1
            rr = rr(1);
        end
        if check == 0
            break
        else
            check = 0; %resetto il check
        end
    end
    
    [Tz_new] = pivoting( Tz,rr,cc ) ;
    
    T_new = Tz_new ;
    
    if gg == 2
         Tz_1 = Tz_new ;
    end


end

