function [alfa,FI,LAMBDA]=alfamaxcontrol(ny,INB,TA,alfamax,FI,LAMBDA,...
    nupt,xj,OUTB)
% 
%
% Control on maximum load factor alfa<alfamax
%
%
    modent=OUTB(ny+1);
    ny2=ny+2;
    if nupt>0
        for i=1:ny
            if INB(i)==0
                alfaold=TA(i,1);
                indexalfa=i;
            end
        end
    end
%
% Computation of alfa at the end of the present loading step
%
%
          if nupt==0
                alfa=xj;
          end
 
          if nupt>0
                step=TA(indexalfa,ny2)*xj;
                alfa=alfaold+TA(indexalfa,ny2)*xj;
          end
 %   
  if alfa>alfamax
      if nupt>0
        beta=(alfamax-TA(indexalfa,1))/TA(indexalfa,ny2);
        alfa =alfamax;
        IN_BASIS_INDECES=INB;
        SOLMAX=TA(1:ny,1)+beta*TA(1:ny,ny2);
      end
      if nupt==0
          alfa=alfamax;
          SOLMAX=TA(1:ny,1)+alfa*TA(1:ny,ny2);
      end
%
%  
     for i=1:ny
         inb=INB(i);
         if inb<=ny
             if inb>0
                FI(inb)=SOLMAX(i);
             end
         end
         if inb>ny
             LAMBDA(inb-ny)=SOLMAX(i);
         end
     end
%
       if nupt>0
            if modent<=ny
                FI(modent)=-beta*TA(ny,ny2);
            end
            if modent>ny
                LAMBDA(modent-ny)=-beta*TA(ny,ny2);
            end
       end
%       
   end %if alfa   
%
end