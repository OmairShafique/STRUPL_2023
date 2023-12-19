function [T_new, Err, T_newp] = PLCP_lastmod(T, last_mod)
%Risolve un problema parametrico per trovare i modi che scaricano
%garantendo l'attivazione dell'ultimo modo attivato

Err = 0 ;

%Effettuo la prima trasformazione pivotale facendo entrare in base il
%parametro
cc = find(T(1,:)==last_mod/10000) ;
mm= T(:,2)./T(:,cc) ;
mm(1)= 100000000 ;  %riga di intestazione
mm(end) = [] ;

if last_mod < 0
    ll = find(T(:,1)==-last_mod) ;
    mm(ll) = 100000000 ; ;
end

p = max(mm(mm<0)) ;
if isempty(p) == 1
    display 'NO solution PLCP - Errore'
    Err = 1 ;
    T_new = T ;
end
if Err == 0
    rr = find(mm==p) ;  %riga elemento Pivot
    if size(rr,1)>1
        rr = rr(1);
    end
    
    %Se un modo LambdaB vuole uscire dalla base,seleziono una diversa riga pivot.
    check = 0 ;  %il check mi serve per interrompere il ciclo quando necessario
    while T(rr,1)< 0 && T(rr,1)>-10000  
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
    
%     %Se l'ultimo modo è di tipo PhiB, devo impedire che la rispettiva
%     %LambdaA esca di base.
%     if last_mod < 0
%         if T(rr,1) == -last_mod/10000
%             mm(rr) = 100000000 ;
%             p = max(mm(mm<0)) ;
%             rr = find(mm==p) ;  %riga elemento Pivot
%             if size(rr,1)>1
%                 rr = rr(1);
%             end
%         end
%     end

    [T_new] = pivoting( T,rr,cc ) ;
    T_newp = T_new;
    
    rr_p = rr ;
end


if Err == 0
    %Proseguo con le altre trasformazioni fino a che la Fi coniugata del
    %parametro esce dalla base
    for yy = 1 : 10000
        T = T_new ;
        if T(1,cc) == last_mod
            for y = 1 : size(T,1)
                if T(y,1) == last_mod/10000
                    plcp_sol = 1 ;
                    break
                end
            end
            if plcp_sol == 1
                break
            end
        elseif T(1,cc)>=10000 || T(1,cc)<=-10000
            cc = find(T(1,:)==T(1,cc)/10000) ;
        else
            cc = find(T(1,:)==T(1,cc)*10000) ;
        end
        mm= T(:,2)./T(:,cc) ;
        mm(1)= 10000000 ;  %riga di intestazione
        
        if last_mod < 0
            ll = find(T(:,1)==-last_mod) ;
            mm(ll) = 100000000 ; ;
        end
        
        p = max(mm(mm<0)) ;
        if isempty(p) == 1
            display 'Unbounded solution - Errore'
            Err = 1 ;
            T_new = T ;
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
            if T(j,2) == 0 && T(j,cc)<0
                rr = j ;
                break
            end
        end
        
%         %Se l'ultimo modo è di tipo PhiB, devo impedire che la rispettiva
%         %LambdaA esca di base.
%         if last_mod < 0
%             if T(rr,1) == -last_mod/10000
%                 -last_mod/10000
%                 mm(rr) = 100000000 ;
%                 p = max(mm(mm<0)) ;
%                 rr = find(mm==p) ;  %riga elemento Pivot
%                 if size(rr,1)>1
%                     rr = rr(1);
%                 end
%                 break
%             end
%         end
        
        %Se un modo LambdaB vuole uscire dalla base,seleziono una diversa riga pivot.
        check = 0 ;  %il check mi serve per interrompere il ciclo quando necessario
        while T(rr,1)< 0 && T(rr,1)>-10000  
            check = 1 ;
            mm(rr)=0 ; %passo a selezionare una nuova riga
            p = max(mm(mm<0)) ;
            if isempty(p) == 1
                display 'Errore rapporto > 0 Lambda B - unbounded solution'
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
        
        if rr == rr_p
            display 'Parametro fuori base - Errore'
            Err = 1 ;
            T_new = T ;
            break
        end
        [T_new] = pivoting( T,rr,cc ) ;
    end
end
end



