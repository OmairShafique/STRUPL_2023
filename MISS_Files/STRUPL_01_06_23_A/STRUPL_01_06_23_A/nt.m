
function [NT,ny,HS,nno,NOT,theta,FI]=nt(ne,ragnetto, ...
        icrk,NCRNPL,ncrl,twonotch)
%    
        theta=0;
%
% ragnetto example
%
%
 if ragnetto==1
   ny=5;
   nno=2
%   nno=0;
   ny=ny+nno;
% 
% ny= total number of plastic modes, including the no tension-compression
% modes;
% nno= total number of no tension-compression modes
% total number of standard modes = ny-nno
% ny=6;total number of plastic modes;
% nno=1 number of notension-compression modes;
% matrix NOT(nyno,2)= matrix which associates standard modes with the
% corresponding no tension-compression mode, FOR RAGNETTO 
% NOT(1,1)=4 (compression softening of bar 3)
% NOT(1,2)=5(no compression mode after softening of mode 4)
%
% Initialization of matrix HS(ny,ny) of hardening or softening
%
    for i=1:ny
        for j=1:ny
            HS(i,j)=0;
        end
    end
%
% Inizialization of matrix NOT
%
        if nno==0
            NOT(1,1)=0
        end
            if nno>0
                for i=1:nno
                 for j=1:2
                    NOT(i,j)=0;
                 end
                end
            end
            
% Definition of matrix NT and matrix HS
%
%
% nrow= number of mode less or equal to ny 
%
%
        nrow=0;
        ncol=0;
        for i=1:ne
%
% Bar 1
%
            if i==1
                %
                % mode in tension
                %
                nrow=nrow+1;
                NT(nrow,i)=1;
            %    
            % mode in compression
            %
                nrow=nrow+1;
                NT(nrow,i)=-1;
                        
            end
%
% Bar 3
%
            if i==3
            %
            % mode in tension
            %
                nrow=nrow+1;
                NT(nrow,i)=1;
            %
            % mode in compression
            %
                nrow=nrow+1;
                NT(nrow,i)=-1;
            %    HS(nrow,nrow)=0;
                HS(nrow,nrow)=-0.2;
            %
            % mode in compression with zero resistance
            %
              nrow=nrow+1;
              NT(nrow,i)=-1;
              NOT(1,1)=4;
              NOT(1,2)=5;
           %
            end
%
% Bar 6
%
            if i==6
           %
           % mode in tension
           %
            nrow=nrow+1;
            NT(nrow,i)=1;
           % HS(nrow,nrow)=0;
            HS(nrow,nrow)=-0.2;
           %
           % mode in tension with zero resistance
           %
              nrow=nrow+1;
              NT(nrow,i)=1;
              NOT(2,1)=6;
              NOT(2,2)=7;
           % 
            end
        end % end for i on elements
    end % end if ragnetto 
%
% cracking problem: 
% icrk=1;
% ncrl= number of cracking lines;
% NCRNPL(ncrl,1)= number of cracking nodes per cracking line
% NOSESPERCRL(ncrl,max(NCRNPL(i),i=1:ncrl))
% modes no tension are activated only in corrispondence with the cut off 
% mode;
% Number of crecking modes = 4x8=32    
% 
  if icrk==1
        irow=0;
        jcol=0;
        inot=0;
        nno=0;
        NOT=0;
        coulomb=1;
        notension=0;
      if twonotch==1
            coulomb=0;
            notension=1
      end
      if coulomb==1
%
% Coulomb friction angle as radians
%
% definition of matrix N'(ny,tncrn*ndf) for Coulomb + tension cut off 
% according to normality rule in the local reference system + hardening
% matrix;
%
 theta =0.523598;
  for i=1:ncrl
     for ii=1: NCRNPL(i)
       irow=irow+1;
       jcol=jcol+1;
  %     inot=inot+1;
  % mode 1 cut off 
       NT(irow,jcol)=1;
       NT(irow,jcol+1)=0;
  % mode 1 twin
  %     NT(irow+1,jcol)=1;
  %     NT(irow+1,jcol+1)=0;
  %
  %       NOT(inot,1)=irow;
  %       NOT(inot,2)=irow+1;
  %
  % mode 2 
  %
       NT(irow+1,jcol)=sin(theta);
       NT(irow+1,jcol+1)=cos(theta);
  %     NT(irow+2,jcol)=sin(theta);
  %     NT(irow+2,jcol+1)=cos(theta);
  % mode 2 twin
  %     NT(irow+3,jcol)=sin(theta);
  %     NT(irow+3,jcol+1)=cos(theta);
  %
  %     NOT(inot+1,1)=irow+2;
  %     NOT(inot+1,2)=irow+3;
  % mode 3
       NT(irow+2,jcol)=sin(theta);
       NT(irow+2,jcol+1)=-cos(theta);
 %      NT(irow+4,jcol)=sin(theta);
 %      NT(irow+4,jcol+1)=-cos(theta);
  % mode 3 twin
  %     NT(irow+5,jcol)=sin(theta);
  %     NT(irow+5,jcol+1)=-cos(theta);
  %
  %     NOT(inot+2,1)=irow+4;
  %     NOT(inot+2,2)=irow+5;
  %
  % mode 1 
      HS(irow,irow)=0;
  %    HS(irow,irow)=-1000;
  % mode 1 twin
  %     HS(irow+1,irow+1)=0;
  % mode 2
        HS(irow+1,irow+1)=0;
  %      HS(irow+2,irow+2)=-1000;
  % mode 2 twin
  %     HS(irow+3,irow+3)=0;
  % mode 3
        HS(irow+2,irow+2)=0;
  %      HS(irow+4,irow+4)=-1000;
  % mode 3 twin
  %      HS(irow+5,irow+5)=0;
  %
  %      inot=inot+2;
        jcol=jcol+1;
  %      irow=irow+5;
        irow=irow+2;
      end %for ii
    end  % for i
     ny=irow;
  %   nno=ny/2;
      end %if coulomb
%
% End Coulomb material and start no tension only
%
  if notension==1
      irow=0;
      jcol=0;
      inot=0;
    for i=1:ncrl
     for ii=1: NCRNPL(i)
       irow=irow+1;
       jcol=jcol+1;
       inot=inot+1;
%       
% mode 1: onnly tension 
%
       NT(irow,jcol)=1;
       NT(irow,jcol+1)=0;
% mode 1 twin
       NT(irow+1,jcol)=1;
       NT(irow+1,jcol+1)=0;
%
       NOT(inot,1)=irow;
       NOT(inot,2)=irow+1;
  %
  %
  % mode 1 
  %     HS(irow,irow)=0;
      HS(irow,irow)=-40476;
      if irow==1
          HS(1,1)=-20238;
      end
      
      if irow==41
          HS(41,41)=-20238;
      end
  % mode 1 twin
       HS(irow+1,irow+1)=0;
  %
          jcol=jcol+1;
          irow=irow+1;
      end %for ii
    end  % for i
     ny=irow;
     nno=ny/2;
   end %if notension
  
 %
 % end no tension material
 %
  end %if icrack
  for j=1:ny
      FI(j)=0;
  end 
%
end

    