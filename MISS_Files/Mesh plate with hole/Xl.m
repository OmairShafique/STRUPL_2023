function xyl = Xl(s,Domain)

global O P1 P2 P3 P4 P5 CMP ;

switch Domain
    case 1
     
        x = P1(1)+(P2(1)-P1(1))*s ;
        y = P1(2)+(P2(2)-P1(2))*s ;
        
    case 2
        
        x = P5(1)+(P4(1)-P5(1))*s ;
        y = P5(2)+(P4(2)-P5(2))*s ;    
        
end

xyl = [x ; y] ;