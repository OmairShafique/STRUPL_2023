function [TA,piv,nupt,YDOT]=pivotran(TA,ny,INB,nupt,toldeta,YDOT, alfamax);
%
% Initial scalars
%
% display ('in pivotran')
%
% Number of Pivot Transformation nupt
%
% display ('Number of Pivot Transformation')
%
nupt=nupt+1
%
ny2=ny+2;
ny1=ny+1;
nym1=ny-1;
%
% Pivot element
%
piv=TA(ny,ny2);
%
TA(ny,ny2)=1/piv;
%
% Elements of the pivot column
%
pivcol=TA(1:nym1,ny2);
TA(1:nym1,ny2)=pivcol*1/piv;
%
% Elements of the NY row
%
pivrow=TA(ny,1:ny1);
TA(ny,1:ny1)=-pivrow*1/piv;
%
% All elements of the matrix TA, a part last column and last row
%
TA(1:nym1,1:ny1)=(TA(1:nym1,1:ny1)-TA(1:nym1,ny2)*pivrow);
display ('out of pivotran')
end