function [FI,LAMBDA,INB,OUTB,alfa,nupt,tol,YDOT]=initialization(ny,BI,R)
%
disp('Sono entrato nella funzione initialization:')
%
%
% display ('in initialization')
tol=.001
alfa=0
nupt=0
FI=R
LAMBDA(ny,1)=0
YDOT(ny,1)=0
INB(ny,1)=0;
OUTB(1,ny+1)=0;
%
%
%
INB(1,1)=1;
OUTB(1,1)=0;
OUTB(1,2)=ny+1;
for i=2:ny
    INB(i)=INB(i-1)+1;
    OUTB(1,i+1)=OUTB(1,i)+1;
end

display ('out of initialization')
end
  