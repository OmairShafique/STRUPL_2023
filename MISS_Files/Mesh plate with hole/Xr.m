function xyr = Xr(s,Domain)
global O P1 P2 P3 P4 P5 CMP ;

switch Domain
    
    case 1
      
        x = CMP(1)+(P3(1)-CMP(1))*s ;
        y = CMP(2)+(P3(2)-CMP(2))*s ;
        
    case 2   

        x = CMP(1)+(P3(1)-CMP(1))*s ;
        y = CMP(2)+(P3(2)-CMP(2))*s ;
          
end

xyr = [x ; y] ;