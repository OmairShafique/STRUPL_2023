function [nupt,TA,pivotindex,piv]=pivotran(TA,ny,nupt,maxasii,tol,flag3)
%,
% display ('in pivotran')
%
%
        pivotindex=0;
        nupt=nupt+1;
%
        ny2=ny+2;
        ny1=ny+1;
        nym1=ny-1;
%
% Pivot element
%
        piv=TA(ny,ny2);
      if abs (piv)<tol %correzione del 28/11-05/12
            pivotindex=1; %correzione del 28/11
            display ('stop in pivotran: pivot=0')%correzione del 28/11
            return %correzione del 28/11
      end %correzione del 28/11
%    
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
%
% All elements of the matrix TA, a part last column and last row
%
        TA(1:nym1,1:ny1)=(TA(1:nym1,1:ny1)-TA(1:nym1,ny2)*pivrow);
%
%               
%display ('out of pivotran')
end