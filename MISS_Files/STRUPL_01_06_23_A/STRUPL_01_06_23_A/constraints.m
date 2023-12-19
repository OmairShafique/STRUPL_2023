
function[tndof,tnfdof,tncdof,SC,TF,TC,TFM,UC]=constraints(nn,ndf,ncn,...
    ragnetto,fourtriangles,twonotch)
%
% display ('enter in constraints')
% It generates matrix TF(tnuncdof, tncondof) where
% tndof=total number of uncontrained structure dof
% tncdof= total number of constrained dof
% tnfdof= total number of free dof
% ncn=number of constrained nodes
%
%      
% Ragnetto problem:
% Total number of unconstrained structure dof=12
% Total number of constraint structure dof =5
% Total number of free dof= 12-5=7
% Total number of master dof=4
% total number of slave nodes =3
%   
%
    for i=1:ncn
       for j=1:15
            SC(i,j)=0;
       end
    end
       if ragnetto==1
%
% Ragnetto example
% dof 1,3,5,9,11 are fixed therefore the corresponding lines are alla zero
% i.e. do not require any intervention;
% dof 2,4,6 are slave of master dof 1,2: t e coefficient woth respwct
% column 1 is 1 while with respect 2 is given by the distance of the slave
% dof with the dof 2, assuming positive 2 if anticlockwise;
% a nalogous reasoning is repeated for slave dof 8,9,10 with respect 
% master 3 and 4.
% 
% Matrx SC(ncn,15) (Structural Constraints) of data concerning constrained 
% nodes are defined;direction X1 coincides with the vertical direction.
% Col 1 = number of the constrained node;
% Col 2-7 node dof (=0 free; =1 constrained)
% Col 8-9 angles of the local reference system;
% Col 10-15 prescribed displacements in the local reference system
%
            SC(1,1)=1;
            SC(1,2)=1;
            SC(2,1)=3;
            SC(2,2)=1;
            SC(3,1)=5;
            SC(3,2)=1;
            SC(4,1)=9;
            SC(4,2)=1;
            SC(5,1)=11;
            SC(5,2)=1;
%
% vectors FDOF( tnfdod,1) (numbering of the free dof) and CDOF (tncdof,1)
% (numbering of the constraned dof) are defined;
%
    [CDOF,FDOF,tndof,tnfdof,tncdof,UC]=consdof(nn,ndf,ncn,SC);
%
% Topological matrix TC(tndof,tncdof) of the constraned dof and TF
% (tndof,tnfdof) of the free dof are determined;
%
    [TC,TF]=tmcs(tndof,tncdof,tnfdof,CDOF,FDOF);

    tnfsdof=3;
    tnfmdof=tnfdof-tnfsdof;
%
% if tnfsdof=0 it means that there are no slave dof and therefore TFM=I
%
    
         if tnfsdof>0
%               
% Defintion of the master and slave nodes dof among the free dof. For the
% ragnetto example:
% 1. u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12 are the unconstrained dof;
% 2. u1,u3,u5,u9,u11 are the constraned dof;
% 3. u2,u4,u6,u7,u8,u10,u12 are the free dof;
% 4. u2,u4,u8,u10 are the master dof;
% 5. u6,u7,u12 are the slave dof;
%
% Equations:
%
% Nodal free displacements UF(tnfdof)=TFM(tnfdof,tnfmdof)*UFM
% Nodal reactions  R(tncdof)=(TC)'*KS*TF*TFM*UM+(TC)'*KS*TC*UC
% Nodal master dis. UFM=((TFM)'*KFF*TFM)-1((TFM)'*FF-(TFM)'*TF*KS*TC*UC)
%
% definition of topological matrix TFM which connects the free dof to the
% master dof:
                      TFM(1,1)=1;
                      TFM(2,1)=1;
                      TFM(3,1)=1;
                      TFM(4,1)=1;
                      TFM(5,1)=0;
                      TFM(6,1)=0;
                      TFM(7,1)=0;
                      TFM(1,2)=3;
                      TFM(2,2)=1;
                      TFM(3,2)=-1;
                      TFM(4,2)=0;
                      TFM(5,2)=0;
                      TFM(6,2)=0;
                      TFM(7,2)=0;
                      TFM(1,3)=0;
                      TFM(2,3)=0;
                      TFM(3,3)=0;
                      TFM(4,3)=0;
                      TFM(5,3)=1;
                      TFM(6,3)=1;
                      TFM(7,3)=1;
                      TFM(1,4)=0;
                      TFM(2,4)=0;
                      TFM(3,4)=0;
                      TFM(4,4)=0;
                      TFM(5,4)=3;
                      TFM(6,4)=1;
                      TFM(7,4)=-1;
         end % if tnfsdof
      end % if ragnetto
%
% Four triangles example
%
    if fourtriangles==1
        SC(1,1)=1;
        SC(1,2)=1;
        SC(1,3)=1;
        SC(2,1)=3;
        SC(2,2)=0;
        SC(2,3)=1;
    end
%
% Two notches example
%
    if twonotch==1
        SC(1,1)=1;
        SC(1,2)=0;
        SC(1,3)=1;
        SC(2,1)=2;
        SC(2,2)=0;
        SC(2,3)=1;
        SC(3,1)=3;
        SC(3,2)=0;
        SC(3,3)=1;
        SC(4,1)=4;
        SC(4,2)=1;
        SC(4,3)=1;
        SC(5,1)=5;
        SC(5,2)=0;
        SC(5,3)=1;
        SC(6,1)=6;
        SC(6,2)=0;
        SC(6,3)=1;
        SC(7,1)=7;
        SC(7,2)=0;
        SC(7,3)=1;
    end %if twonotch==1
%    
%   
%
    if ragnetto==0
        [CDOF,FDOF,tndof,tnfdof,tncdof,UC]=consdof(nn,ndf,ncn,SC);
%
        [TC,TF]=tmcs(tndof,tncdof,tnfdof,CDOF,FDOF);
%
%
%
       tnfsdof=0
          if tnfsdof==0
            for i=1:tnfdof
                for j=1: tnfdof
                    TFM(i,j)=0;
                    if i==j
                    TFM(i,j)=1;
                    end
                end
            end
          end
    end
   end












