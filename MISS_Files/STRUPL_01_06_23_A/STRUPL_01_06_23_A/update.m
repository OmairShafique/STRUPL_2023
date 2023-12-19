function [alfa,INB,OUTB,FI,LAMBDA,flag3,FLAG1,LASTPAIR,TA,flag2]=...
    update(INB,OUTB,ny,TA,nno,NOT,FLAG1,indexi,lambdas,flag3,LASTPAIR,flag2)

%
                
%
%
        if flag3==0
            if flag2==1
                TA(ny,1)=TA(ny,1)+lambdas;
                flag2=0;
            end
        end
%       
%
% display ('entering in updateindeces')
%
        ny1=ny+1;
%
%--------------------------------------------------------------------------
%-------------------------------------------------------------------------- 
%------ Updating INB and OUTB vectors--------------------------------------
%
    imode=INB(ny,1);
    jmode=OUTB(1,ny+1);
    INB(ny,1)=jmode;
    OUTB(1,ny+1)=imode;
% 
%----- Calculation of alfa and vectors FI and LAMBDA----------------------
%
%
        for i=1:ny
            FI(i,1)=0;
            LAMBDA(i,1)=0;
        end
        alfa=0;
        for i=1:ny
            kk= INB(i);
            if kk==0
                alfa=TA(i,1);
            end
            if kk > ny
                LAMBDA(kk-ny,1)=TA(i,1);
            end
            if kk>0
             if kk<=ny
                 FI(kk,1)=TA(i,1);
             end
            end
         end
%
%--------------------------------------------------------------------------
%-------Updating lambda value of the no tension mode which substitutes-----
%-------the softening mode-------------------------------------------------
%
% The lamda of a no tension mode which has substituded a softening mode is
% updated with the value of the lambda of the preceeding softening lambda
% value;

%                lambdas=lambdas
%                inbny= INB(ny)
%           for i=1:nno
%               ntwmlam=NOT(i,2)+ny
%                if inbny==ntwmlam
%                    TA(ny,1)=TA(ny,1)
%               TA(ny,1)=TA(ny,1)+lambdas
%                end
%            end




%-----------------update LASTPAIR, flag3,flag2 and FLAG1-----------------------------------
% 
             lastoutmode=OUTB(ny+1)
%
% it is the last mode which has been activated in the last pivot 
% transformation: the pair(one Fi and one Lambda) present contemporanealy 
% in OUTB are memorized in a vector named LASTPAIR(1,2). 
% The first component represent the last activated mode while the second 
% is the mode already out of basis.
%
%--------------------------------------------------------------------------
           if lastoutmode<=ny
               if flag3==0
                LASTPAIR(1,1)=lastoutmode;
                LASTPAIR(1,2)=lastoutmode+ny;
               end
             for i=1:nno;
                if lastoutmode==NOT(i,1)
                   FLAG1(i,2)=0;
                end
                if lastoutmode==NOT(i,2)
                   for j=1:ny+1
                      outmod=OUTB(j);
                        if outmod==lastoutmode-1
                           flag3=1;
                           FLAG1(i,1)=1;
                          end % if outmod
                    end % for j
                end % if lastoutmode=NOT(i,2)
              end % for i
           end % if lastoutmode<=ny
      if lastoutmode>ny
           if flag3==0
                LASTPAIR(1,1)=lastoutmode;
                LASTPAIR(1,2)=lastoutmode-ny;
% ------------------------------------------------------------------------                
                indfi=lastoutmode-ny;
                 for m=1:nno
                    if NOT(m,1)==indfi 
                     FLAG1(m,2)=1;
                    end
                 end
           end
%-------------------------------------------------------------------------
           if flag3==1
               flag2=1
               flag3=0
               for i=1:nno
                   if (lastoutmode-ny)==NOT(i,1)
                    FLAG1(i,1)=1
                   end %if
               end %for i
           end % if flag3
      end%if lastoutmode > ny
%
%            
% -------------------------------------------------------------------------
% display ('Out from upindeces')       
end  
