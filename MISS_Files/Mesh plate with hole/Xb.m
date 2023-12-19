function xyb = Xb(s,Domain)
global R ;
r = R ;
global O P1 P2 P3 P4 P5 CMP ;

switch Domain
    case 1
       
        x = O(1)+r*cos(pi/4*s) ;
        y = O(2)+r*sin(pi/4*s) ;
           
    case 2

        x = O(1)+r*cos(pi/4*s) ;
        y = O(2)+r*sin(pi/4*s) ;
       
end

xyb = [x ; y] ;

