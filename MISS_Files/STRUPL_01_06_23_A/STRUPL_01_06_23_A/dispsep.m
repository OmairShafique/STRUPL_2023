function [UMEP,UGAMMAP,CRDIS]= dispsep (KSFMM1,FSFM,alfa,KS,TS,TF,...
    TC,UC,TFM,LAMBDA,NT,DS,BS,icrk,KE,TCR,RCR,ndf,ncrl,NCRNPL,NODESPERCRL,nn)
%
%
%   computes the free master elastic plastic displacements UMEP due to the 
%   free master nodal loads FSFM and given prescribed nodal prescribed 
%   displacements UC and palstic deformations LAMBDA in the case of plastic
%   behaviour (icrk=0) or cracking behavior (icrk=1). In the cracking
%   behavior UMEP collects the nodal displacement on the  Gamma- side of the
%   cracks, which do not move for cracking, while UGAMMAP the displacements 
%   of the Gamma+ side, which moves.
%
%  
%
%    
    UGAMMAP=0;
    CRDIS=0;
%
%   
     if icrk==0
        % UMEP1=TF*TFM*KSFMM1*FSFM*alfa
        % UMEP2=TF*TFM*KSFMM1*TFM'*TF'*KS*TC*UC
        % UMEP3=TF*TFM*KSFMM1*TFM'*TF'*TS'*BS'*DS*NT'*LAMBDA
        % UMEP4=TC*UC
        % UMEP=UMEP1-UMEP2+UMEP3+UMEP4
        % UMEP=TF*TFM*KSFMM1*(FSFM*alfa-TFM'*TF'*KS*TC*UC...
        % +TFM'*TF'*TS'*BS'*DS*NT'*LAMBDA)+TC*UC;
        % UMEP=TF*TFM*KSFMM1*(FSFM*alfa-TFM'*TF'*KS*TC*UC...
%         +TFM'*TF'*TS'*BS'*DS*NT'*LAMBDA)+TC*UC;
     
      UMEP=KSFMM1*(FSFM*alfa-TFM'*TF'*KS*TC*UC...
         +TFM'*TF'*TS'*BS'*DS*NT'*LAMBDA);
     
     end

    if icrk==1    
        UMEP=TF*TFM*KSFMM1*(FSFM*alfa-TFM'*TF'*KS*TC*UC...
         +TFM'*TF'*TS'*KE*TCR'*RCR'*NT'*LAMBDA)+TC*UC;
%     

         [UGAMMAP,CRDIS]=discrack(UMEP,ndf,LAMBDA,NT,RCR,ncrl,NCRNPL,...
                            NODESPERCRL,nn);
    end
%
%
end