function xyt = Xt(s,Domain)

global O P1 P2 P3 P4 P5 CMP ;

switch Domain
    case 1   

        x = P2(1)+(P3(1)-P2(1))*s ;
        y = P2(2)+(P3(2)-P2(2))*s ;
        
    case 2

        x = P3(1)+(P4(1)-P3(1))*s ;
        y = P3(2)+(P4(2)-P3(2))*s ;    

end

xyt = [x ; y] ;
