

function [LAMBDA,FI,FLAG1]=lcpy(ny,A,RS,BI,alfamax,KSFMM1,FSFM,TF,KS,TS,...
         TC,UC,TFM,NT,DS,BS,nno,NOT,maxasii,icrk,KE,TCR,RCR,ndf,ncrl,...
         NCRNPL,NODESPERCRL,nn,FS,inslc,FI,FLAG1,tol,twonotch,ragnetto)
%
% 
% The parametric procedure(parameter alfa is the positive scalar which 
% amplifies proportionally all load in the loading step delta t)describes 
% the finite step elastic plastic problem.  
% It is conisidered as given data the following vectors and matrix:scalar
% ny as the number of activable plastic modes, vector FIZERO(ny) (vector of the
% ny initial plastic potentials),matrix A(ny,ny), BI(ny) vector of the 
% elastic stresses due to the load increments during the step pojected 
% on the normal of the activable plastic modes 
% dei carichi nel passo proiettati sui modi plastici BI, 
%
%               Parametric Linear Complementarity Problem
%
%             FI(t+delta t)= FI(t)+ A*delta Lambda-alfa*BI 
%           FI >=0; delta Lambda>=0; (FI)t(delta Lambda)=0;
%
% Mechanical Remarks: 
% the above formulation will be able to solve:
% (i) the elastic-plastic-softening static holonomic problem both for a
% unique proportional loading condition and a sequence of proportional 
% loading stages, like for ciclic loadings. 
% (ii) the elastic plastic softening dinamic problem in which,
% at every time step,a proportional loading condition is solved.
% First implementation will deal with the one step proportional static
% loading.
% Numerical aspects to be underlined:
% (iii) Generally speaking, the first step of the PLCP coincides with 
% the solution of an LCP as:
%   FIDOT=A*LAMBDADOT-BIDOT; 
%   FIDOT>=0; LAMBDADOT>=0;(FIDOT)t*LAMBDADOT=0
% because vector FI(t) has some component equal to zero; 
% (iV) on the contrary, at the first step of the first proportional loading 
%      stage, vector FI(t)=RS (vector of plastic reisstances) which can be
%      assumed all positive and therefore the PLCP can have a trivial first
%      step solution performing a simple pivot transformation by entering 
%      into the basis the load factor alfa end exiting from the basis the 
%      first componebt of the vector FI which becomes equal to zero.
% PLCP can be described via a TABLEAU form trougth matrix TA which connects
% , at the starting procedure, variables FI (variables in the Basis) with
% parameter alfa and variables LAMBDA.
% The rules are:
% 1. at each transformation (step), one non basic variable enters into the
%   Basis and one variable in the Basis goes out of the Basis;
% 2. the parameter alfa increases until the problem is definite
%   positive;in special cases alfa remains constant during a Tabeau
%   transformation (multiplicity of solutions at constant load);
%   in the case of undefinite problem ALFA decreases;
% 3. all variables are always positive or zero;
% 4.2 complementary (dual) variables satisfy always the ortogonality
%   condition FI(i)*LAMBDA(i)=0;
% Matrix TA(ny,ny+2) is made of ny rows and ny+2 columns. 
% First column contains the values of the in Basis variables; initially 
% it coincides with the initial values of vector FI=RS >= 0.
% Second column it represents, initially, vector -BI corresponding to
% parameter alfa out of Basis.
% All other columns coincide with the colums of matrix A.
% Basis variables are ny while Out of Basis variables are ny+1. 
% Parameter alfa will enter into the Basis at first pivotal transformation
% and it will remain in Basis for the entire loading procedure.
% The loading procedure terminates when the prescribed value of alfa 
% has been reached or because a limt load has been attained. In the case of
% softening behavior the procedure will stop when the laod has reach a
% local maximum; the procedure can restart to follow the descending path
% untill until the load multiplier ALFA has attained the prescribed
% decreased value or the problem has returned to be definite positive.
% Vectors LAMBDA and FI are updated at each pivotal transformation,
% together with the load factor alfa.
% At each pivotal transformation one non basic varialble (=0) enters into
% the Basis (>0) and one basic variable (>0) goes out of the Basis
% (decoming =0).
% For the moment the assumption of "non degenenancy" is adopted, i.e. it is
% supposed that only one dual variables LAMBDA(i),FI(i) can be
% contemporaneausly =0.
% We know theat "degenerancy" of the solution will be able to deal with
% multiplicity of solutions, i.e. with "equilibrium bifurcation points".
% At the first pivot transformation, variable alfa enters into the Basis
% and a variable FI exit from the Basis; at the second transformation, we
% find out of the Basis one cople of dua vasiables, i.e. LAMBDA(i) and
% FI(i); therefore at the second transformation the LAMBDA(i) variable is
% going to enter into the Basis.
% Therefore the procedure assumes always as a general rule that the
% variable that is goint to enter into the Basis is the dual variable which
% has just exit from the Basis.
% In the Scientific Literature this rule is called "Restricted Basis rule".
% A part from particular situation( when the matrix A of the active modes
% become singular or undefined)the generic pivot transformation is
% represented by the following occurence:
% (i) one active mode (FI(i)=0 and LAMBDA(i)>0) become unactive, then first
% LAMBDA(i)becomes =0 and then FI(i) becomes>0, or
% (ii) one mode mode i unactive(i.e. FI(i)>0 and LAMBDA(i)=0) becomes
% active then first FI(i) becomes =0 and then LAMBDA(i) becomes >0.
% From the numerical point of view, in order to keep track of the variables
% which are in the Basis and the variables which are out of the Basis at
% each Pivot Transformation, the fowwing vectors are created and updated at
% each pivot transformation:
% Column INB(i), with i=1,ny which list the code of the basic variables.
% The adopted code index is the following;
% code =0 corresponds to parameter alfa;
% code 1 to ny list the variables FI(i);
% code ny+1 to 2ny list the variables LAMBDA8i).
% The codes of the non basic variables are listed in the row vector
% OUTB(1,1+ny).
% Data are made available in Workspace memory: number of plastic modes ny, 
% matrix A, vectors RS, BI;
%
% Print input values for ny,A,BI,RS, alfamax,nupt
% Inizialization of unknown vectors FI,LAMBDA,column vector INB (indices of 
% basic bariables) and row vector OUTB (indices of out of basis variables)
% Initialization of vector column INB and vector column OUTB
%
%--------------------------------------------------------------------------
%--------------------INITIALIZATION OF VARIABLES---------------------------
   disp('.....Entering into initialization  from lcpy function.....')
%--------------------------------------------------------------------------
%
    [FI,LAMBDA,INB,OUTB,alfa,nupt,LAMBDADOT,indexjstandard,FLAG1,...
       LASTPAIR,flag2,lambdas]=initialization(ny,BI,RS,nno,inslc,FI,FLAG1);
%
%--------------------------------------------------------------------------    
%---------------INITIALIZATION OF TABLEAU TA------------------------------- 
% 
        [TA,indexj]=tableauzero(FI,BI,A,tol,ny);
        
%        COL1=TA(:,1)
%
%--------------------------------------------------------------------------
   disp('.....Terminated initialization  from lcpy function.....')
%--------------------------------------------------------------------------
%---------------------------NOTATION---------------------------------------
%
% nupt= number of pivot pransformations
% indexi : row index is the row number of matri TB which is going to exit
% from the Basis;
% INB(i)= number of unknown variable (FI if <=ny or LAMBDA if>ny) going out
% from the Basis;
% indexj= column index is the column number of matrix TB which is going to
% enter into the Basis;
% OUTB(j)= number of the variable (FI if <=ny or LAMBDA if >ny) entering
% into the Basis.
% It starts a sequence of functions, i.e
% function outbasis,
% function rearrange,
% function pivotran,
% function updateindeces,
% function inbasis
% which perform pivot transformations until the final load factor is
% reached or before if a collapse situation is determined.
% "function outbasis" computes row i (indexi)
% of matrix TB which is going to exit from the Basis. 
%--------------------------------------------------------------------------
% 
%
   nusomo=0;
   flag3=0;
     while alfa<alfamax
       Step_Number="Step Number="+(nupt+1)
%
%   failure mechanism
%
%---------------------FAILURE MECHANISM------------------------------------        
%        if indexfail==1
%--------------------------------------------------------------------------
 disp('.....Entering into failuremech function.....')
%--------------------------------------------------------------------------
%
           [icollapse]=failuremech(ny,TA,indexj,INB,OUTB,alfa)
%        
%
%--------------------------------------------------------------------------
 disp('.....Terminated failuremech function.....') 
%--------------------------------------------------------------------------
            if icollapse==1
                 return
            end 
%
%--------------------------------------------------------------------------
 disp('.....Entering into outbasis function.....')
%--------------------------------------------------------------------------
%
%--------COMPUTES THE MODE NUMBER EXITING FROM THE BASIS-------------------
%
   [xj,indexi,lambdas,TA,INB,OUTB,imode,jmode,indexfail,indstop]=outbasis(TA,...
       indexj,ny,nno,NOT,INB,OUTB,maxasii,tol,flag3,nusomo,FLAG1,flag2,...
       lambdas,nupt);
   if indstop==1
       return
   end

%--------------------------------------------------------------------------
 disp('.....Terminated outbasis function.....')
%--------------------------------------------------------------------------
%
                     imode=imode
%
%--------------------------------------------------------------------------
%------------REARRANGE MATRIX TA AND VECTORS INB AND OUTB------------------
%
%   Vectors INB,OUTB and matrx TA are rearranged in order to change row ì
%   with row ny and Column j with column ny+2 (last column of the Tableau TA.
%
%--------------------------------------------------------------------------
  disp('.....Entering into rearrange function.....')
%--------------------------------------------------------------------------
%
        [TA, INB,OUTB]=rearrange(TA,INB,OUTB, indexi,indexj,ny);
%
%--------------------------------------------------------------------------
  disp('.....Terminated rearrange function.....')
%--------------------------------------------------------------------------
%
%----------------CONTROL ON THE MAXIMUM LOAD FACTOR ALFA-------------------
% Contol on alfa<alfamax: if alfa greather then alfamax, then the solution
% for alfa =alfamax is computed and the program stops.
%
%--------------------------------------------------------------------------
    disp('.....Entering into alfamaxcontrol function.....')
%--------------------------------------------------------------------------

        [alfa,FI,LAMBDA]=alfamaxcontrol(ny,INB,TA,alfamax,FI,LAMBDA,...
            nupt,xj,OUTB);
        if alfa==alfamax
            return
        end
%        
%-------------------------------------------------------------------------- 
   disp('.....Terminated alfamaxcontrol function.....')
%--------------------------------------------------------------------------
 
%--------------------------------------------------------------------------
%-----------------PIVOT TRANSFORMATION-------------------------------------
% Pivot transformation on TA matrix with pivot the element TA(ny,ny+2), it
% modifies the ny row of matrix TA(ny+1 first terms), it modifies the ny+2
% column of matrix TA (First NT-1 terms, and the remaining terms of marix TA
% (ny-1 rows and ny+1 columns) with matrix multiplications.
% If the Pivot term is =0 the program send to the failure mechanism function.
%
%--------------------------------------------------------------------------
    disp('.....Entering into pivotran function.....')
%--------------------------------------------------------------------------
%
       [nupt,TA,pivotindex,piv]=pivotran(TA,ny,nupt,maxasii,tol,flag3);
%                            pivot=piv
%   
%--------------------------------------------------------------------------
    disp('.....Terminated pivotran function.....')
%--------------------------------------------------------------------------   
          if pivotindex==1 %correzione del 29/11
        % [LAMBDADOT,FI]=failuremech(ny,TA,indexj,INB,OUTB)
           return % correzione del 29/11
          end % if pivotindex correzione del 29/11
%          
% vectors INB and OUTB after the pivot transformation: INB(ny) is
% the row index which exit from the Basis and OUTB(ny+2) is the column index
% which enters into the Basis. INB is the column vectors which lists the
% indeces of FI( from 1 to ny) and OUTB lists the LAMBDA (from ny+1 to 2ny)
% vectors.Number 0 identifies the load parameter which stays always into 
% the Basis, a part the first step when it goes from out of Basis into the
% Basis.
%
%--------------------UPDATING----------------------------------------------
       disp('.....Entering into update function.....')
%--------------------------------------------------------------------------
%   
 [alfa,INB,OUTB,FI,LAMBDA,flag3,FLAG1,LASTPAIR,TA,flag2]=...
  update(INB,OUTB,ny,TA,nno,NOT,FLAG1,indexi,lambdas,flag3,LASTPAIR,flag2);
                            alfa=alfa
                            INBASIS=INB'
                            OUTBASIS=OUTB
                           % NOT=NOT
                            FLAG1=FLAG1
                            fag3=flag3
                            flag2=flag2
                           LAMBDAt=LAMBDA'
                           FIt=FI'
                           INBASEVECTORt=TA(:,1)'
%
%
%--------------------------------------------------------------------------
    disp('.....Terminated update function.....')
%--------------------------------------------------------------------------
%--------------------Control on "creasy" loop------------------------------
                    if nupt==15
 disp('--------------------Max number of steps----------------------')
                       return
                    end
%--------------------------------------------------------------------------
%
%------COMPUTATION OF THE MODE ENTERING IN THA BASIS----------------------- 
% Choice of the new column j entering into the Basis
% The choice of the variable out of Basis which is going to enter into the 
% the Basis is determined by the following rules:
% First step: variable alfa, second column of matrix TA;
% generic step: dual variable of the variable which has exit from the Basis
% If the variable which is going to exit from the Basis,say row i, has 
% a code index INB(i)<= ny (i.e it is a variable FI, then we have to 
% find the dual variable in row vector OUTB with  code index 
% OUTB(j)= INB(i)+ny while if OUTB(i) >ny (i.e. it is a LAMBDA variable
% we have to find j such that OUTB(j)=INB(i)-ny.
%
%--------------------------------------------------------------------------
     disp('.....Entering into inbasis function.....')
%-------------------------------------------------------------------------- 
                inrow=OUTB(1,ny+1);
%
%  
        [indexj,colj,modej]=inbasis(OUTB,ny,inrow,nno,NOT,...
                             INB,flag3,indexjstandard,TA,LASTPAIR);
%
%--------------------------------------------------------------------------
     disp('.....Terminated inbasis function.....')
%--------------------------------------------------------------------------
            if indexj==0
                return
            end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------        
% Computes the displacements due to the external loads, prescribed nodal
% displacements and "plastic strains" lambda; 
%
%--------------------------------------------------------------------------
%----------------------DISPLACEMENTS COMPUTATION---------------------------
    disp('.....Entering into dispsep function.....')
%--------------------------------------------------------------------------
%
                        if flag2==0
%                            
% if flag2=1 it means that a void transformation has been just performend
% by setting =0 a lambda variable correspong to the softening mode which is
% no more active without the activation of the corresponding no tension
% mode a nd therefore with a lambda vector with no mechanical sense.
%
% 
%--------------------------------------------------------------------------
      [UMEP,UGAMMAP,CRDIS]=dispsep(KSFMM1,FSFM,alfa,KS,TS,TF,...
 TC,UC,TFM,LAMBDA,NT,DS,BS,icrk,KE,TCR,RCR,ndf,ncrl,NCRNPL,NODESPERCRL,...
 nn);
%   
    if twonotch==1
            ucontrolR=UGAMMAP(528*2+2)
            ucontrolM=UGAMMAP(525*2+2)
            ucontrolL=UGAMMAP(522*2+2)
    end
    if ragnetto ==1
        ucontrol=UMEP(3)
    end
%                        
%--------------------------------------------------------------------------
    disp('.....Terminated dispsep function.....')
%--------------------------------------------------------------------------
%
%-------------------COMPUTATION OF REACTIONS-------------------------------
% Computes reactions in elastic plastic range
%
%--------------------------------------------------------------------------
    disp('.....Entering into reactionsep function.....')
%--------------------------------------------------------------------------
%
        [RSEP]=reactionsep(UMEP,UC,TC,TF,KS,TFM,TS,BS,DS,NT,LAMBDA,...
            icrk,KE,TCR,RCR,FS);
%
%--------------------------------------------------------------------------
    disp('.....Terminated reactionsep function.....') 
%
                        end %if flag2=0
%--------------------------------------------------------------------------
%
% The program goes back to the computation of the indexi of the Basis element 
% which is going to exit, i.e. to function outbasis   disp ('alfa=alfamax')
%   
%--------------------------------------------------------------------------
  end % while

 end
