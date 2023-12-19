% 
% The parametric procedure(parameter ALFA is the positive scalar which 
% amplifies proportionally all load in the loading step delta t)describes 
% the finite step elastic plastic problem.  
% It is conisidered as given data the following vectors and matrix:scalar NY
% as the number of activable plastic modes, vector FIZERO(NY) (vector of the
% NY initial plastic potentials),matrix A(NY,NY), BI(NY) vector of the 
% elastic stresses due to the load increments during the step pojected 
% on the normal of the activable plastic modes 
% dei carichi nel passo proiettati sui modi plastici BI, 
%
%               Parametric Linear Complementarity Problem
%
%             FI(t+delta t)= FI(t)+ A*delta Lambda-ALFA*BI 
%           FI >=0; delta Lambda>=0; (FI)t(delta Lambda)=0;
%
% Mechanical Remarks: 
%
% the above formulation will be able to solve:
% (i) the elastic-plastic-softening static holonomic problem both for a
% unique proportional loading condition and a sequence of proportional 
% loading stages, like for ciclic loadings. 
% (ii) the elastic plastic softening dinamic problem in which,
% at every time step,a proportional loading condition is solved.
% 
% First implementation will deal with the one step proportional static
% loading.
%
% Numerical aspects to be underlined:
% 
% (iii) Generally speaking, the first step of the PLCP coincides with 
% the solution of an LCP as:
%
%   FIdot=A*LAMBDAdot-BIdot; 
%   FIdot>=0; LAMBDAdot>=0;(FIdot)t*LAMBDAdot=0
%
% because vector FI(t) has some component equal to zero; 
%
% (iV) on the contrary, at the first step of the first proportional loading 
%      stage, vector FI(t)=R (vector of plastic reisstances) which can be
%      assumed all positive and therefore the PLCP can have a trivial first
%      step solution performing a simple pivot transformation by entering 
%      into the basis the load factor ALFA end exiting from the basis the 
%      first componebt of the vector FI which becomes equal to zero.
%
% PLCP can be described via a TABLEAU form trougth matrix TA which connects
% , at the starting procedure, variables FI (variables in the Basis) with
% parameter ALFA and variables LAMBDA.
%
% The rules are:
% 1. at each transformation (step), one non basic variable enters into the
%   Basis and one variable in the Basis goes out of the Basis;
% 2. the parameter ALFA increases until the problem is definite
%   positive;in special cases ALFA remains constant during a Tabeau
%   transformation (multiplicity of solutions at constant load);
%   in the case of undefinite problem ALFA decreases;
% 3. all variables are always positive or zero;
% 4.2 complementary (dual) variables satisfy always the ortogonality
%   condition FI(i)*LAMBDA(i)=0;
%
% Matrix TA(NY,NY+2) is made of NY rows and NY+2 columns. 
% First column contains the values of the in Basis variables; initially 
% it coincides with the initial values of vector FI=R >= 0.
% Second column it represents, initially, vector -BI corresponding to
% parameter ALFA out of Basis.
% All other columns coincide with the colums of matrix A.
% Basis variables are NY while Out of Basis variables are NY+1. 
%
% Parameter ALFA will enter into the Basis at first pivotal transformation
% and it will remain in Basis for the entire loading procedure.
% The loading procedure terminates when the prescribed value of ALFA 
% has been reached or because a limt load has been attained. In the case of
% softening behavior the procedure will stop when the laod has reach a
% local maximum; the procedure can restart to follow the descending path
% untill until the load multiplier ALFA has attained the prescribed
% decreased value or the problem has returned to be definite positive.
%
% Vectors LAMBDA and FI are updated at each pivotal transformation,
% together with the load factor ALFA.
% At each pivotal transformation one non basic varialble (=0) enters into
% the Basis (>0) and one basic variable (>0) goes out of the Basis
% (decoming =0).
% For the moment the assumption of "non degenenancy" is adopted, i.e. it is
% supposed that only one dual variables LAMBDA(i),FI(i) can be
% contemporaneausly =0.
% We know theat "degenerancy" of the solution will be abple to deal with
% multiplicity of solutions, i.e. with "equilibrium bifurcation points".
% At the first pivot transformation, variable ALFA enters into the Basis
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
%
% From the numerical point of view, in order tokeep track of the variables
% which are in the Basis and the variables which are out of the Basis at
% each Pivot Transformation, the fowwing vectors are created and updated at
% each pivot transformation:
% Column INB(i), with i=1,NY which list the code of the basic variables.
% The adopted code index is the following;
% code =0 corresponds to parameter ALFA;
% code 1 to NY list the variables FI(i);
% code NY+1 to 2NY list the variables LAMBDA8i).
% The codes of the non basic variables are listed in the row vector
% OUTB(1,1+NY).
%
% 
% Data are made available in Workspace memory: number of plastic modes NY, 
% matrix A, vectors R, BI;
%
     load ragnettodata.mat
%
% Print input values for NY,A,BI,R, ALFAMAX,NUPT
%
% Inizialization of unknown vectors FI,LAMBDA,column vector INB (indices of 
% basic bariables) and row vector OUTB (indices of out of basis variables)
%
% 
%
% Initialization of vector column INB and vector column OUTB
% 
%
    [FI,LAMBDA,INB,OUTB,ALFA,NUPT,TOL,YDOT]=initialization(NY,BI,R)
%
% Initialization Tableau TA 
% 
    [TA,INDEXJ,TOLDETA]=tableauzero(R,BI,A,TOL,NY);
%
% Notation:
% NUPT= number of pivot pransformations
% Index i : row index is the row number of matri TB which is going to exit
% from the Basis;
% INB(i)= number of unknown variable (FI if <=NY or LAMBDA if>NY) going out
% from the Basis;
% Index j= column index is the column number of matrix TB which is going to
% enter into the Basis;
% OUTB(j)= number of the variable (FI if <=NY or LAMBDA if >NY) entering
% into the Basis.
%
% It starts a sequence of functions, i.e
% function outbasis,
% function rearrange,
% function pivotran,
% function updateindeces,
% function inbasis
% which perform pivot transformations until the final load factor is
% reached or before if a collapse situation is determined.
%
%
% "function outbasis" computes row i (INDEXI)
% of matrix TB which is going to exit from the Basis. 
%
%
% Loop on the undefinite pivotal trnsformations
%
   while ALFA<ALFAMAX
   [DELTA,INDEXI]=outbasis(TA,INDEXJ,NY);
%
% Vectors INB,OUTB and matrx TA are rearranged in order to change row ì
% with row NY and Column j with column NY+2 (last column of the Tableau TA.
%
    [TA, INB,OUTB]=rearrange(TA,INB,OUTB, INDEXI,INDEXJ,NY);
%
% control on failure mechanism
%
    [YDOT]=failuremech(NY,TA,INDEXJ);
    if FLAG ==1
    return
    end 
%
%
%
% Pivot transformation on TA matrix with pivot the element TA(NY,NY+2), it
% modifies the NY row of matrix TA(NY+1 first terms), it modifies the NY+2
% column of matrix TA (First NT-1 terms, and the remaining terms of marix TA
% (Ny-1 rows and NY+1 columns) with matrix multiplications.
% If the Pivot term is =0 the program stops with a message of error.
%
    [TA,PIV,NUPT,YDOT]=pivotran(TA,NY,INB,NUPT,TOLDETA,YDOT)
%
% Update of vectors INB and OUTB after the pivot transformation: INB(NY) is
% the row index which exit from the Basis and OUTB(NY+2) is the colun index
% which enters into the Basis. INB is the column vectors which lists the
% indeces of FI( from 1 to NY) and OUTB lists the LAMBDA (from NY+1 to 2NY)
% vectors.Number 0 identifies the load parameter which stays always into 
% the Basis, a part the first step when it goes from out of Basis into the
% Basis.
%
    [INB,OUTB,inrow]=updateindeces(INB,OUTB,NY);
%

  %
  % Choice of the new column J entering into the Basis
  % The choice of the variable out of Basis which is going to enter into the 
  % the Basis is determined by the following rules:
  % First step: variable ALFA, second column of matrix TA;
  % generic step: dual variable of the variable which has exit from the Basis
  % If the variable which is going to exit from the Basis,say row i, has 
  % has a code index INB(i)<= NY (i.e it is a variable FI, then we have to 
  % find the dual variable in row vector OUTB with  code index 
  % OUTB(j)= INB(i)+NY while if OUTB(i) >NY (i.e. it is a LAMBDA variable
  % we have to find j such that OUTB(j)=INB(i)-NY.
  %
  % update the value of ALFA
  %
        [ALFA]=upalfa(TA,INB,NY);
  %
  %
        [INDEXJ]=inbasis(OUTB,NY,inrow);
  %
 %
 % The program goes back to the computation of the INDEXI of the Basis element 
 % which is going to exit, i.e. to function outbasis      disp ('ALFA=ALFAMAX')
    end
  

