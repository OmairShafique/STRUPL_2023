

function [xj,indexi,lambdas,TA,INB,OUTB,imode,jmode,indexfail,indstop]=outbasis(TA,indexj,ny,nno,...
    NOT,INB,OUTB,maxasii,tol,flag3,nusomo,FLAG1,flag2,lambdas,nupt)
%
% Computes the index i (variable indexi row of matrix TA) of the 
% basic element which is going to exit out of the Basis;
% if all terms of column j are positive then a situation of unbounded
% solution is determined and the procedure will stop with the associated
% plastic mechanism
%
% display ('entering in outbasis')
%    
%
%-----------Case of flag3=1 index void transformation----------------------
%     
% the softening mode lambda is in the Basis, 
% The lambda softening mode (modei) goes out of the basis
% with constant alfa (void transformation) entering in basis 
% the correspong FI of the softening mode set =0.
%
             imode=0;
             indexfail=0;
             indstop=0;
%             
% Computation of mode number=modej corresponding to indexj     
%
               jmode=OUTB(indexj-1);
%        
%           
    if flag3==1 
        xj=0;
%
    modei=OUTB(ny+1)-1+ny;
            for m=1:ny
                if INB(m)==modei
                    indexi=m;
                    lambdas=TA(indexi,1);
                    TA(indexi,1)=0;
                end %if INB
            end % for m
 end %if flag3=1
%
%
%-------------   End of the case flag3=1-----------------------------------
%
%-----------Computation of flag1: index to avoid mode FI control(=1)-------
  if flag3==0
%
%--------flag2=1 means to add lambdas to no tension lambda mode------------
    if flag2==0
        lambdas=0;
    end
%--------------------------------------------------------------------------%
      xj=10^20;
      indexi=0;
    for i=1:ny
      flag1=0;
      iv=INB(i,1);
        if iv<=ny 
           if nno>0
            for j=1:nno
                inds=NOT(j,1);
                indno=NOT(j,2);
                    if iv==inds
                        flag1=FLAG1(j,1);
                    end %if inds
                    if iv==indno
                        flag1=FLAG1(j,2);
                    end %if indno
            end %for j
          end %if nno
        end %if iv<=ny
%
% Control on FI and Lambda positivness
%
                if flag1==0 
                   % if TA(i,indexj)<-1E-05
                        rap(i)=-TA(i,1)/TA(i,indexj);
                         if rap(i)>0
                            if rap(i)<xj 
                                xj=rap(i);
                                indexi=i;
                                imode=INB(i);
                            end %if rap<xj
                         end %if rap>0
                   % end %if TA
               end %if flag1
    end %for i
 end %if flag3=0
%
% Control on positivness >=0 of the j column corresponding to the mode
% entering in the Basis
%
                if nupt>0
%
%
% Search for the row index ialfa corresping to variable alfa in the Basis
                    ialfa=0;
                    for k=1:ny
                        if INB(k)==0
                            ialfa=k;
                        end
                    end
                    if ialfa==0
              display ('alfa=0 out of the basis')
              indstop=1
                    end
              return
%                    
% Control on the zero value of the matrix term T(ialfa,indexj): if it is
% zero it means multiplicity of solution i.e. the configuration corresponds
% to a mechanism (not always plastic i.e. satisfying normality rule);
%
               if TA(ialfa,indexj)<tol
                    if TA(ialfa,indexj)>-tol
                        indexfail=1;
                    end 
                end
%                        
% control on the positivness of the the remaining indexj column
%
                        if indexfail==1
                            for m=1:ny
                                if TA(m,indexj)<-tol
                                    indexfail=0;
                                end
                            end
                        end
           end %if nupt
%                                
% display ('exiting from outbasis')
%
 end % end of the function




  