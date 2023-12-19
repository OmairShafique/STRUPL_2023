function[ncrl,tncrn,NCRNPL,NODESPERCRL,BETA,RCR,INFLEN]=crlines(fourtriangles,COORD)
%
% defines ncrl= number of cracking lines;
% defines vector NCRNPL of the Number of Nodes Per Cracking Line;
% defines tncrn= total number of cracking nodes;
%
% defines NODESPERCRL, i.e. the nodes per each cracking line;
% computes vector BETA of the angles between the cracking line and axis x1;
% comptes matrix RCRS(2*tncn,2*tncn) of the rotation matrix from axis x1 
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
                % calaculation of the influence lenght of the nodes os
                % cracling line 1
                %
                % starting node
                
                
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
            % computes the total number of nodes along various cracking
            % lines;
            %
                 tncrn=0;
                for i=1:ncrl;
                 tncrn=tncrn+NCRNPL(i);
                end
            %    ny=3*tncrn;
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
               
          end %if fourtringles
        
    end