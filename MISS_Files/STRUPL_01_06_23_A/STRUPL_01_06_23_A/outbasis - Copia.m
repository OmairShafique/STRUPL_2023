function [delta,indexi]=outbasis(TA,indexj,ny)
%
% Computes the index i (variable INDEXI) of the basic element which is 
% going to exit out of the Basis;
% if all terms of column j are positive then a situation of unbounded
% solution is determined and the procedure will stop with the associated
% plastic mechanism
%
% display ('entering in outbasis')
    
    indexi=0
    delta=10E23;
    for i=1:ny
        if TA(i,indexj)<0
            
             rap(i)=-TA(i,1)/TA(i,indexj);
            
               if rap(i)<delta 
                    delta=rap(i);
                    indexi=i;
               end
         end
     end
% display ('exiting from outbasis')
end
    