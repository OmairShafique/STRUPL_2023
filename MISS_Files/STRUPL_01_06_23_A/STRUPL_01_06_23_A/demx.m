function [DE]= demx(nel,nste,elt,ragnetto)
   if ragnetto==1
% ragnetto case
     for i=1:nste
      for j=1:nste
        DE(i,j)=0;
          if i==j
            DE(i,j)=1; 
            
                if nel==5
                DE(i,j)=0.5 ;
                end
                if nel==6
                DE(i,j)=0.5;
                end
           end
          end
     end
   end
    if elt==3
% 2D  plane stress problem
          E=18000;
          nu=0.2;
      
      
        for i=1:nste
            for j=1:nste
                DE(i,j)=0;
            end
        end
        A=E/(1-nu^2);
        DE(1,1)=1;
        DE(2,2)=1;
        DE(1,2)=nu;
        DE(2,1)=DE(1,2);
        DE(3,3)=0.5*(1-nu);
        DE(:,:)=DE(:,:)*A;
%    
    end
end