function [TCR] = topcrack(ne,ndf,ncrl,NCRNPL,NODESPERCRL,ELNODES,ONESELPERCRN,...
    ncolones,tncrn)
% Input:
% ndf=number of node dof;
% ncrl= number of potential cracking lines;
% NCRNPL= number of cracking nodes per line;
% NODESPERCRL= node numbers per cracking line;
%
% Output:
% TCR(tncrn*ndf,nex3xndf) topology matrix (0-1) which transforms element nodal
% displacements into cracking node displacements;
%
    tnrow=tncrn*2;
    tncol=6*ne;
        for i=1:tnrow
            for j=1:tncol
                TCR(i,j)=0;
            end
        end
      icrnode=0;
%
% irow is the row index
%
         for i=1:ncrl
          nnopcrl=NCRNPL(i);
%
            for j=1:nnopcrl
               icrnode=icrnode+1;
%                
                  for q=1:ncolones
                    iel=ONESELPERCRN(icrnode,q);
%
% if che controlla l'uscita dal loop sul numero di elemento da una parte
% della fessura;
%
                      if iel>0
%
% Loop sul numero di nodi dell'elemento;
%
                        for n=1:3
                         if ELNODES(iel,n)==NODESPERCRL(i,j) 
                             irow=(icrnode-1)*ndf;
                             jcol=(iel-1)*6+(n-1)*ndf;
                                for h=1:ndf;
                                  for k=1:ndf
                                    if k==h
                                       TCR(irow+h,jcol+k)=1;
                                    end
                                  end
                                end
                         end %for n

                      end % if iel
                  end % for q
               
              end %for j
        end %for i
% 
   end