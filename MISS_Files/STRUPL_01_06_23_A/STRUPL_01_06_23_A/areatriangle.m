function [areatr]=areatriangle(nel,ELNODES,COORD)
%
% Computes the area of a triangle with the Erone formula
%
    na=ELNODES(nel,1);
    nb=ELNODES(nel,2);
    nc=ELNODES(nel,3);
%
    XA(1)=COORD(na,1);
    XA(2)=COORD(na,2);
    XB(1)=COORD(nb,1);
    XB(2)=COORD(nb,2);
    XC(1)=COORD(nc,1);
    XC(2)=COORD(nc,2);
%
   ab1=abs(XA(1)-XB(1));
   ab2=abs(XA(2)-XB(2));
   ab=sqrt(ab1^2+ab2^2);
%    
    ac1=abs(XA(1)-XC(1));
    ac2=abs(XA(2)-XC(2));
    ac=sqrt(ac1^2+ac2^2);
%
    bc1=abs(XB(1)-XC(1));
    bc2=abs(XB(2)-XC(2));
    bc=sqrt(bc1^2+bc2^2);
% 
%   
    p=(ab+ac+bc)/2;
    areatr=sqrt(p*(p-ab)*(p-ac)*(p-bc));
%
end

 
