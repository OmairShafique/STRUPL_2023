
function [BE,ndfe]=bemx(nel,nste,ndf,elt,COORD,ELNODES);
    if elt==0
% ragnetto problem 
        ndfe=2;
        for i=1:nste
            for j=1:ndfe
                if j==1
                    BE(i,j)=-1;
                else
                    BE(i,j)=1;
                end
            end
        end
    end
%
    if elt==3
%
%   plane stress triangular plane stress 3 nodes element  
        ndfe=ndf*3;
%
% Inizializie matrix BE, matrix B of the element
        for i=1:3
            for j=1:6
                BE(i,j)=0;
            end
        end

%
%  Genrate matrix X(6,6), function of the nodal coordinates
%  
        for i=1:6
            for j=1:6
                X(i,j)=0;
            end
        end
%
% computes matrix X
%
        nn1=ELNODES(nel,1);
        nn2=ELNODES(nel,2);
        nn3=ELNODES(nel,3);
        X(1,1)=1;
        X(1,2)=COORD(nn1,1);
        X(1,3)=COORD(nn1,2);
        X(2,4)=1;
        X(2,5)=COORD(nn1,1);
        X(2,6)=COORD(nn1,2);
        X(3,1)=1;
        X(3,2)=COORD(nn2,1);
        X(3,3)=COORD(nn2,2);
        X(4,4)=1;
        X(4,5)=COORD(nn2,1);
        X(4,6)=COORD(nn2,2);
        X(5,1)=1;
        X(5,2)=COORD(nn3,1);
        X(5,3)=COORD(nn3,2);
        X(6,4)=1;
        X(6,5)=COORD(nn3,1);
        X(6,6)=COORD(nn3,2);
%
% Compute the inverse of matrix X
%
        XM1=inv(X);
%
% Compute matrix B for the element nel
        BE(1,:)=XM1(2,:);
        BE(2,:)=XM1(6,:);
        BE(3,:)=0.5*(XM1(3,:)+XM1(5,:));
  end

        
 
end

