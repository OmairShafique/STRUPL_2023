
function [icollapse]=failuremech(ny,TA,indexj,INB,OUTB,alfa);
%
% display ('entering in failuremech')
% collapse mechanism
% 
% failure mech verifies if the configuration (matrix TA) is at collapse by
% controlling 2 conditions:
% Condition N.1: scalar TA(i,indexj) with i corresponding to load factor
% alfa should be =0;
% Condition N.2: all other terms of column TA(i,indexj) should be >= 0;
% if both conditions are satisfied the procedure stops otherwisw it conues
% with function outbasis.m.
%
% Inintializatiion
% 
            icollapse=0;
            ialfa=0
%
% verification of the collapse configuration
%
            for i=1:ny
                if INB(i)==0
                    ialfa=i;
                end
            end
            if ialfa>0
              if abs(TA(ialfa,indexj))<=1E-04
                icollapse=1;
                for j=1:ny
                    if icollapse==1
                        if TA(j,indexj)<0
                            if abs(TA(j,indexj))>1E-04
                                icollapse=0;
                            end %if abs
                        end %if TA
                    end %if icollapse
                end %for j
              end %if abs
            end %if ialfa
%
            if icollapse==1
%
% Initialization
%
                for i=1:ny
                    FI(i,1)=0;
                    LAMBDADOT(i,1)=0;
                end
                    iout=OUTB(indexj-1)
                if iout<=ny
                    TACJ(:,1)=TA(:,indexj)
                end
                if iout>ny
                    LAMBDADOT(iout-ny,1)=1;
                end
                for i=1:ny
                    in=INB(i);
                        if in<=ny
                            if in>0
                                FI(in,1)=TA(i,1);
                            end %if in>0
                            if in==0
                                alfa=TA(i,1);% aggiunta del 29-11
                            end %if in=0
                        end %if in<=ny
                        if in>ny
                            indexydot=in;
                            LAMBDADOT(in-ny,1)=TA(i,indexj);
                        end %if in>ny
                end %for i
                alfa=alfa
                LAMBDADOT=LAMBDADOT
                FI=FI
            end %if collapse
  %          
  %display ('exit fron failuremech')
  %
   end