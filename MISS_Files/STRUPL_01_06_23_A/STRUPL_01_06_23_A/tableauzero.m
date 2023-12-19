function [TA,indexj]=tableauzero(FI,BI,A,tol,ny)
 display (' in tableauzero')
%
%
%
[TA]=[FI,-BI,A];
%
%
%  Compute the tolerance on the zero value of the diagonal term and the 
%  pivotal column of A-1 matrix. This tolerance is used to estimate a zero
%  value of the matrix A-1.
%  The adopted rule is: TOLam1=10-04xmin[ABS(1/Aii)
%  toldeta=(TA(1,3))*tol
%
%
%        
 indexj=2;
 display (' out from tableauzero')
end

