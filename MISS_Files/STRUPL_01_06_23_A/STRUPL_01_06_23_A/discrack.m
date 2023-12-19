
function [UGAMMAP,CRDIS]=discrack(UMEP,ndf,LAMBDA,NT,RCR,ncrl,NCRNPL,NODESPERCRL,nn)
%
% Computes the vector UGAMMAP whch lists the structure nodal displacements
% in the global reference system due to cracking with respct the negative
% gamma lines, i.e. the line which remain fixed during the cracking phase.
% If a node is not iterested to cracking the UGAMMAP conponents are zero.
%
% Computation of the cracking nodes displacements relative to the global
% reference system:
%
    ndftot=nn*ndf;
    for i=1:ndftot
        UGAMMAP(i,1)=0;
    end
%
    CRDIS=RCR'*NT'*LAMBDA;
%
    icrnode=0;
    for i=1:ncrl
      for j=1:NCRNPL(i)
          icrnode=icrnode +1;
          nodenum=NODESPERCRL(i,j);
          iugammap=(nodenum-1)*ndf;
          icrdis=(icrnode-1)*ndf;
          for k=1:ndf
             UGAMMAP(iugammap+k,1)=CRDIS(icrdis+k,1);
          end %for k
      end % for j
     end %for i
%
    for i=1:ndftot
        UGAMMAP(i,1)=UGAMMAP(i,1)+UMEP(i,1);
    end
%      deltapl=UMEP(528*2+2)
%      deltatot=UGAMMAP(528*2+2)
end