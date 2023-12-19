  function [delta,Reaction]=solve(KK,fg,bc)
% SOLVE --> Solve static FE-equations considering boundary conditions.
%-------------------------------------------------------------
% INPUT:   K  : global stiffness matrix, dim(K)= nd x nd, nd : number of dof's
%          F  : global load vector, dim(f)= nd x 1
%          bc : boundary condition matrix
%               dim(bc)= nbc x 2, nbc = number of b.c.'s
%               bc = (gdl pd), pd = prescribed displacemnt
%
% OUTPUT:  u : displacement solution including boundary values
%          R : reaction force vector
%              dim(u)=dim(R)= nd x 1 
%-------------------------------------------------------------
  if nargin==2 ; 
     delta=KK\fg;  
  elseif nargin==3;
     [nnd,nnd]=size(KK);
     fdof=[1:nnd]';  %free dof
     
     delta=zeros(size(fdof));
     
     pdof=bc(:,1);  %prescribed dof
     pd=bc(:,2);    %prescribed displacement
     fdof(pdof)=[]; 

     s=KK(fdof,fdof)\(fg(fdof)-KK(fdof,pdof)*pd);
    
     delta(pdof)=pd;
     delta(fdof)=s;
  end  
     Reaction=KK*delta-fg;
%--------------------------end--------------------------------
