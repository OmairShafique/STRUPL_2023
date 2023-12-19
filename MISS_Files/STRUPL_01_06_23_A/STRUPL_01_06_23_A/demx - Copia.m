function [DE]= demx(nel,nste,elt)
if elt>0
    return
end

   for i=1:nste
      for j=1:nste
        DE(i,j)=0;
        
           if i==j
            DE(i,j)=1; 
            
                if nel==5
                DE(i,j)=0.5 
                end
                if nel==6
                DE(i,j)=0.5
                end
               
            end
        end
    end
end