function[ncrl,tncrn,NCRNPL,NODESPERCRL,BETA,RCR,INFLEN]=...
    Crack_Model(fourtriangles,COORD,twonotch)
%
% defines ncrl= number of cracking lines;
% defines vector NCRNPL of the Number of Nodes Per Cracking Line;
% defines tncrn= total number of cracking nodes;
%
% defines NODESPERCRL, i.e. the nodes per each cracking line;
% computes vector BETA of the angles between the cracking line and axis x1;
% computes matrix RCRS(2*tncn,2*tncn) of the rotation matrix from axis x1 
% in the local v (v1 and v2) reference system, where v1 coincides with the 
% outward normal to the cracking line;
%
        if fourtriangles==1
            ncrl=3;
%
% definition of the nodes of the cracking lines 
%
% cracking line 1
%
                NCRNPL(1)=2;
                NODESPERCRL(1,1)=1;
                NODESPERCRL(1,2)=3;
%
% global reference system and local reference systems are 
% considered positive anticlockwise;
% definition of angle beta of the normal to the
% crackung line is calculated with respect to axis x1; 
% angle beta is computed as beta =-3.14/2+alfa where alfa is
% the angle of the crack direction with respect axis x1;
%
                L2=COORD(3,2)-COORD(1,2);
                L1=COORD(3,1)-COORD(1,1);
                Le2=L1^2+L2^2;
                den=sqrt(Le2);
                alfa=acos(L1/den);
                BETA(1)=-pi/2+alfa;
%
% cracking line 2;
%  
                NCRNPL(2)=3;
                NODESPERCRL(2,1)=1;
                NODESPERCRL(2,2)=5;
                NODESPERCRL(2,3)=4;
                L2=COORD(4,2)-COORD(5,2);
                L1=COORD(4,1)-COORD(5,1);
                Le2=L1^2+L2^2;
                den=sqrt(Le2);
                alfa=acos(L1/den);
                BETA(2)=-pi/2+alfa;
%
% cracking line 3;
%
                NCRNPL(3)=3;
                NODESPERCRL(3,1)=2;
                NODESPERCRL(3,2)=5;
                NODESPERCRL(3,3)=3;
                L2=COORD(2,2)-COORD(5,2);
                L1=COORD(2,1)-COORD(5,1);
                Le2=L1^2+L2^2;
                den=sqrt(Le2);
                alfa=acos(L1/den);
                BETA(3)=-pi/2+alfa;
            %    
         end %if fourtriangles
 %        
 % 
         if twonotch==1
                 ncrl=1;
 %
 % definition of the nodes of the cracking lines 
 %
 % cracking line 1
 %
                NCRNPL(1)=21;
                NODESPERCRL(1,1)=253;
                NODESPERCRL(1,2)=254;
                NODESPERCRL(1,3)=255;
                NODESPERCRL(1,4)=256;
                NODESPERCRL(1,5)=257;
                NODESPERCRL(1,6)=258;
                NODESPERCRL(1,7)=259;
                NODESPERCRL(1,8)=260;
                NODESPERCRL(1,9)=261;
                NODESPERCRL(1,10)=262;
                NODESPERCRL(1,11)=263;
                NODESPERCRL(1,12)=264;
                NODESPERCRL(1,13)=265;
                NODESPERCRL(1,14)=266;
                NODESPERCRL(1,15)=267;
                NODESPERCRL(1,16)=268;
                NODESPERCRL(1,17)=269;
                NODESPERCRL(1,18)=270;
                NODESPERCRL(1,19)=271;
                NODESPERCRL(1,20)=272;
                NODESPERCRL(1,21)=273;
%                
%               
%
                L2=COORD(253,2)-COORD(273,2);
                L1=COORD(253,1)-COORD(273,1);
                Le2=L1^2+L2^2;
                den=sqrt(Le2);
                alfa=acos(L1/den);
                BETA(1)=-pi/2+alfa;
          end %if twonotch
 %        
 % calaculation of the influence lenght of the nodes 
 %
 %
                 tncrn=0;
                for i=1:ncrl;
                 tncrn=tncrn+NCRNPL(i);
                end
 %
 %   Computation of matrix INFLEN(ncrl,max (NCRNPL(i)i=1,ncrl)
 %
                for i=1:ncrl
                    ncrn=NCRNPL(i);
                    for j=1:ncrn
                         if j==1
                            nodenumm1=NODESPERCRL(i,j);
                            nodenump1=NODESPERCRL(i,j+1); 
                         end %if j
                          if j>1
                            if j<ncrn
                                nodenumm1=NODESPERCRL(i,j-1);
                                nodenump1=NODESPERCRL(i,j+1);
                            end %if j
                            end %if j
                          if j==ncrn
                            nodenumm1=NODESPERCRL(i,j-1);
                            nodenump1=NODESPERCRL(i,j);
                          end %if j
                        dx1=abs(COORD(nodenumm1,1)-COORD(nodenump1,1));
                        dx2=abs(COORD(nodenumm1,2)-COORD(nodenump1,2));
                        INFLEN(i,j)=0.5*sqrt(dx1^2+dx2^2);
                     end % for j
                 
                    end %for i
%
% Computation of matric RCR (tncrn*2,tncrn*2)
% of the rotation matrix from local
% to global reference system for all cracking nodes;
%
                    index=0;
               for i=1:ncrl;
                   B=BETA(i);
                   for j=1:NCRNPL(i)
                        index=index+1;
                        RCR(index,index)=cos(B);
                        RCR(index,index+1)=sin(B);
                        index=index+1;
                        RCR(index,index-1)=-sin(B);
                        RCR(index,index)=cos(B);
                   end %for j
               end %for i
   %
   %
         end