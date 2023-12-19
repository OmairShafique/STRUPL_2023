function [FS]=gravityloads(COORD,ELNODES,ELTYPE,nn,ndf,ne,gama,VOL,DIRG)
%
% Computes nodal loads FS due to element gravity loads gama
%
%
      for i=1:nn
          for j=1:ndf
              nrow=(i-1)*ndf+j;
              FS(nrow,1)=0;
          end %for j
      end %for i
%
       for i=1:ne
         if ELTYPE(i)==3
%
% node numbers of elemnt i (triangular CSS
%
            nna=ELNODES(i,1);
            nnb=ELNODES(i,2);
            nnc=ELNODES(i,3);
%
% Nodal coordinates of nodes nna,nnb,nnc
%
            xa1=COORD(nna,1);
            xa2=COORD(nna,2);
            xb1=COORD(nnb,1);
            xb2=COORD(nnb,2);
            xc1=COORD(nnc,1);
            xc2=COORD(nnc,2);
%
% Coordinates of the centroid
%
            xG1=(1/3)*(xa1+xb1+xc1);
            
            xG2=(1/3)*(xa2+xb2+xc2);
%
% Coordinates with respct the centroid
%
            xGa1=xa1-xG1;
            xGa2=xa2-xG2;
            xGb1=xb1-xG1;
            xGb2=xb2-xG2;
            xGc1=xc1-xG1;
            xGc2=xc2-xG2;
%
% Definitin of matrix XBAR
%
            XBAR(1,1)=1;
            XBAR(2,1)=1;
            XBAR(3,1)=1;
            XBAR(1,2)=xGa1;
            XBAR(2,2)=xGb1;
            XBAR(3,2)=xGc1;
            XBAR(1,3)=xGa2;
            XBAR(2,3)=xGb2;
            XBAR(3,3)=xGc2;
 %
 % Inverse of matrix XBAR
 %
            XBARM1=inv(XBAR);
 %
 %  Computation of vectors C1,C2,C3
 %
            C1(1)=XBARM1(1,1);
            C1(2)=XBARM1(2,1);
            C1(3)=XBARM1(3,1);
             
            C2(1)=XBARM1(1,2);
            C2(2)=XBARM1(2,2);
            C2(3)=XBARM1(3,2);
                
            C3(1)=XBARM1(1,3);
            C3(2)=XBARM1(2,3);
            C3(3)=XBARM1(3,3);
 %
 % Computation of nodal element loads equivalent to gravity forces
 %
            
            FS(2*nna-1,1)=FS(2*nna-1,1)+gama*VOL(i)*DIRG(1)*C1(1);
            FS(2*nna,1)=FS(2*nna,1)+gama*VOL(i)*DIRG(2)*C1(1);
            
            FS(2*nnb-1,1)=FS(2*nnb-1,1)+gama*VOL(i)*DIRG(1)*C2(1);
            FS(2*nnb,1)=FS(2*nnb,1)+gama*VOL(i)*DIRG(2)*C2(1);
            
            FS(2*nnc-1,1)=FS(2*nnc-1,1)+gama*VOL(i)*DIRG(1)*C3(1);
            FS(2*nnc,1)=FS(2*nnc,1)+gama*VOL(i)*DIRG(2)*C3(1);
         end %if ELTYPE
       end %for i
    end
    









   