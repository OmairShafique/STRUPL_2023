
%
%   Script Program strupl2023 (STRUctural PLasticity 2023)
%
%   Convenction: all scalars and indices are written with a lower case while
%   vectors and matrices are written with the upper case;
%
% 
%   matrices:
%   COORD(nn,ncord),ELNODES(ne,maxnnperel),NNODEPEREL(ne,maxnnperel),
%   ELSPERNODE(nn,maxelpernode),
%   TS(tneldof,tndof),TF(tndof,tnfdof),TC(tndof,tncdof),TFM(tnfdof,tnmfdof),
%   BS(tnelstress,tneldof),DS(tnelstress,tnelstress),
%   KS(tndof,tndof),KSF(tnfdof,tnfdof),KSFM(tnmfdof,tnmfdof);
%
%   
%   NT(ny,tnelstress),AS(ny,ny),
% 
%
%   vectors:
%   NNODEPEREL(ne),ELTYPE(ne);NSTRESSPEREL(ne)
%   LAMBDA(ny),FI(ny),R(ny);
%
%
%   scalars:
%   nn,ne,tneldof,prtype,ncord,maxnnperel,maxelpernode,tneldof,
%   tndof,tnfdof,tncdof,tnelstress,ny,ncn,nslc
%
%
%   Function geometry defines scalars and matrices explained at the 
%   beginning of the function; 
%   
%
%        load ('Coordinates.txt')
%        COORD=Coordinates
%        load ('Element_matrix.txt')
%         ELNODES=Element_matrix
%         return
%
    [ncn,ndf,COORD,ELNODES, prtype,nn,ne,maxnnperel,maxelpernode, ...
    ELSPERNODE,NNODEPEREL,tneldof,tndof,ELTYPE,NSTRESSPEREL,...
    ragnetto,th,fourtriangles,icrk,nslc,ALFAMAX,gama,DIRG,ncrl,...
    NODESPERCRL,ONESELPERCRN,NCRNPL,RCR,TCR,INFLEN,tol,twonotch]=geometry;
%
%
%   Function topology computes  matrix TS (Structure Topology matrix) which
%   relates the element nodal  dof (6 bar elements*2 dof per elements=12 for
%   "ragnetto" example)and 2 elements*6 dof per elements = 12 dof to the
%   struture nodal unconstrained displacements = 12 for the "ragnetto"
%   example;
%
%
    [TS]=topology(ndf,nn,ELNODES,ne,maxnnperel,NNODEPEREL,ELTYPE);
%
%
%   Function constraints receives as input matrix SC(ncn,15)(where 
%   ncn = number of constrained nodes) which defines data of
%   constraned nodes and the topological matrices TC(tndof,tncdof)and
%   TF(tndof,tnfdof) which relate the structural unconstrained dof to the
%   contrained and free dof respectively. In other words, these toplogical
%   matrices have the role to rearrange the nodal displacement vector in a
%   first part where the free dof are listed and a second final part where
%   the constrained dof are positioned;
%   
%
    [tndof,tnfdof,tncdof,SC,TF,TC,TFM,UC]=constraints(nn,ndf,ncn,...
    ragnetto,fourtriangles,twonotch);
%
%   
%   Function bsmx computes elemnt matrices BS which relate the generalized
%   element deformations to the elements nodal displacements for the entire
%   structure;
%
%
    [BS]=bsmx(ne,nn,ndf,COORD,NNODEPEREL,ELNODES,NSTRESSPEREL,ELTYPE);
%
%
%   Function dsmx computes matrix DS, relating element stress to strain,
%   times the elemnt volume for the case of element constant stress and
%   strain, for alla elements of the structure in a diagonal form;
%
%
    [DS,VOL]=dsmx(NSTRESSPEREL,ne,ELTYPE,ragnetto,ELNODES,COORD,th);
% 
%
%   Function ksmx computes:
%   KE: diagonaql matrix of the element stiffness matrices;
%   KS: structure stiffness matrix including the constraned dof;
%   KSF: structure stiffness matrix referring only to the free dof;
%   KSFM: structure stiffness matrix referring to master free dof; 
%
%
    [KSFM,KSF,KS,KE]=ksmx(TS,TF,BS,DS,TFM);
%    
%
%   Function ksmxm1 Computes the inverse KSFMM1 of the structural stiffness 
%   matrix KSFM, relative to the master free dof;
%
%
    [KSFMM1]=ksmxm1(KSFM);
%
%
%   Computation of topoloogical matrix TRC(tncrn*ndf,ne*3*ndf) which
%   connects the craks relative displacements to the nodal element
%   displacements;
%
    if icrk == 1 ;
                          % se icrk=1 materiale cracking; se icrk=0
                          % materiale plastico. Per ora si prevede che i
                          % modi siano solo di un tipo,
                          % o cracking o plastici.
%      
%
%   Definition of:
%               1. the number of cracking lines "ncrl";
%               2. the total number of cracking nodes "tncrn";
%               3. total number of cracking modes "ny";
%               4. vector "NCRNPL" which gives the numper of cracking nodes
%                  per cracking line;
%               5. matrix "NODESPERCRL" which lists the node number for
%                  each cracking line;
%               6. vector "BETA" with the angles between the normal 
%                   versors to the cracking lines;
%               7. matrix "RCR", diagonal matrix of the transformation of
%                  the global reference system to the local reference system
%                  for each cracking line and each cracking node; 
%               8. matrix "INFLEN" which gives the influence length of each
%                  node for each cracking line;
%               9. matrix TCR(ndf*tncrn,6*ne): toological matrix 0/1 whch
%                  transforms global element displacements into local nodal
%                  cracking dsplacements;
%               
       [ncrl,tncrn,NCRNPL,NODESPERCRL,BETA,RCR,INFLEN]=...
       crlines(fourtriangles,COORD,twonotch);
%
%   definition of matrix ONESELPERCRN with the elements on one side of the 
%   crack for each cracking line and each node of the cracking line;
%
     [ONESELPERCRN,ncolones]=oneselpercrn(ncrl,NCRNPL,NODESPERCRL,COORD,...
              ELNODES,ELSPERNODE,BETA,maxelpernode);
%
%   definiion of the topological matrix TCR which connects the element
%   displacements in the global reference system with the nodal crack 
%   displacements;
%
%
   [TCR] = topcrack(ne,ndf,ncrl,NCRNPL,NODESPERCRL,ELNODES,ONESELPERCRN,...
          ncolones,tncrn);
%
%
%   Yield conditions are defined throught matrix NT for the "ragnetto"
%   example;it lists the yield conditions for tension and compression of bar
%   1 and 3 and only tension for bar 6; 
%   It can be expressed in general by writing the unit outward normals 
%   to piecewise linear yield conditions;
%   In the case of non linear yield conditions, the yield function is
%   linearized locally with reference to the elastic solution, adopting a
%   backward or forward criterion. The final solution is then checked with
%   reference to the respect of the yield condition violation. If the error is
%   greather then the allawable error, a new linearization is added to the
%   original with respect to the final stress state and a new solution,
%   with more constraints, is determined. The procedure continues until the
%   solutiionis is satisfactory;the above procedure will be implemented in
%   a second time.
%
%   
%
     end  % if icrk
%
%
%   Computation of matrix NT, i.e. N transposed;
%
%
    [NT,ny,HS,nno,NOT,theta,FI]=nt(ne,ragnetto,icrk,NCRNPL,ncrl,twonotch);
    
%
%
%   Computation of matrix AS;
%
%
    [AS,maxasii]=am(NT,DS,BS,TS,TF,KSFMM1,TFM,HS,RCR,TCR,ragnetto,...
        icrk,KE,ny);
%
%
%   Computes element plastic resistances R;
%
%
    [RS]=res(ny,ragnetto,icrk,ncrl,NCRNPL,theta,INFLEN,th,twonotch)
%
%
%-----------LOOP ON SEQUENCE OF PROPORTIONAL LOADING CONDITIONS------------
%------------------CONDITION N. 1 SELF WEIGTH------------------------------
%
%------------INITIALIZATION OF VECTOR LAMBDATOT and FLAG1------------------
%--------------------------------------------------------------------------
%
               FLAG1=0;

               for m=1:nno
                    FLAG1(m,1)=0;
                    FLAG1(m,2)=0;
               end
%
               for k=1:ny
                    LAMBDATOT(k,1)=0;
               end
%
        for i=1:nslc
            inslc=i
display ('---------------Loading Conditiion Number-----------------------')  
                            INSLC=inslc
display ('---------------------------------------------------------------')
%
            alfamax=ALFAMAX(i);
%   Function nodal loads defines  nodal load vector FSFM for the specific 
%   cases of ragnetto and twotriangles;
%
%
    [FSFM,FS]=nodalloads(COORD,ELNODES,ELTYPE,nn,ndf,ne,gama,VOL,DIRG,...
             ragnetto,TFM,TF,fourtriangles,inslc,twonotch);
%           
%
%  Computes the elastic displacement vector US relative to the free dof;
%
%   
    [UME]=dispse(KSFMM1,FSFM,TF,KS,TC,UC,TFM);
%
%
%   Computes the elastic reactions vector RS for the costrained dof;
%
%
    [RSE]=reactionse(UME,UC,TC,TF,KS,TFM,FS);
%
%   Computes vector BI, linear function of external loads
%
    [BI]=bi(ny,NT,DS,BS,TS,TF,TFM,KSFMM1,FSFM,ragnetto,icrk,KE,...
        TCR,RCR);
%
%   Solver "finite" elastic plastic (hardening/softening) problem for
%   alfamax;
%
    [LAMBDA,FI,FLAG1]=lcpy(ny,AS,RS,BI,alfamax,KSFMM1,FSFM,TF,KS,TS,...
        TC,UC,TFM,NT,DS,BS,nno,NOT,maxasii,icrk,...
        KE,TCR,RCR,ndf,ncrl,NCRNPL,NODESPERCRL,nn,FS,inslc,FI,FLAG1,tol,...
        twonotch,ragnetto);
%
%------------------------UPDATING OF LAMBDATOT-----------------------------
%
                LAMBDATOT=LAMBDATOT+LAMBDA;
%
%--------------------------------------------------------------------------
%   end of loading conditions
%--------------------------------------------------------------------------
        end
%--------------------------------------------------------------------------
%
    
    
    