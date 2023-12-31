function [FI,LAMBDA,INB,OUTB,alfa,nupt,LAMBDADOT,indexjstandard,FLAG1,...
    LASTPAIR,flag2,lambdas]=initialization(ny,BI,RS,nno,inslc,FI,FLAG1)
%
% disp('Entering into initialization function:')
%
%
%
        lambdas=0;
        flag2=0;
        indexjstandard=0;
        
        alfa=0;
        nupt=0;
%
        if inslc==1
            FI=RS;
        end
%
        LAMBDA(ny,1)=0;
        LAMBDADOT(ny,1)=0;
        INB(ny,1)=0;
        OUTB(1,ny+1)=0;
        LASTPAIR(1,1)=0;
        LASTPAIR(1,2)=0;
%
%   Initialization of vectors INB and OUTB
%
        INB(1,1)=1;
        OUTB(1,1)=0;
        OUTB(1,2)=ny+1;
    for i=2:ny
        INB(i)=INB(i-1)+1;
        OUTB(1,i+1)=OUTB(1,i)+1;
    end
%
% Initialization of vector FLAG1(nno,2);FLAG1(i,1) refers to softening mode
% while FLAG1(i,2) refers to no tension mode;
%
        if nno>0
            if inslc==1
                for i=1:nno
                    FLAG1(i,1)=0;
                    FLAG1(i,2)=1;
                end %for i
            end %if inslc
        end %if nno
%        
%display ('out of initialization function')
end
  