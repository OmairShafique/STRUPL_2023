


function[ncn,ndf,COORD,ELNODES,prtype,nn,ne,maxnnperel,maxelpernode, ...
    ELSPERNODE,NNODEPEREL,tneldof,tndof,ELTYPE,NSTRESSPEREL,...
    ragnetto,th,fourtriangles,icrk,nslc,ALFAMAX,gama,DIRG,ncrl,...
    NODESPERCRL,ONESELPERCRN,NCRNPL,RCR,TCR,INFLEN,tol,twonotch]=geometry;
%
%               Scalars
%
% prtype= 2 plane truss; 
%       = 3 space truss; 
%        =4 plane frame;
%        =5 space frame;
%        =6 plane stress/strain;
%        =7 space brick;
%        =8 plane slab;
%        =9 space shell
%
% ncoor=2 plane problem; =3 space problem; number of nodal coordinates
%
% nn= number of nodes
% ne= number of elements
% ncn= number of constrained nodes
% th=thickness of the plane stress problem
%
% ndf= number of nodal degrees of freedom n.d.f.
%   ndf=2 for plane truss and plane stress and strain problems ;
%   ndf=3 for plane frames, space brick and plate problems ;
%   ndf=6 for space frame and shell problems
%
% maxnnperel= max number of nodes per elements
% tneldof = total number of element dof
% tndof= total number of unconstrained structure dof
% tnelstress= total number of elemnt stress components
% 
% maxelpernode= Max number of elements around a node
% maxnnperel= Max number of nodes per element
% nslc= number of sequential loading conditions
%
%             Vectors and matrices
% ALFAMAX(nslc) max values of sequential loading conditions
% COORD(nn,ncoor)  nodal coordinates;
% ELNODES(ne, maxnnperel)nodes of each  element;
% ELSPERNODE (nn,maxelpernode) elements around each node;
% NNODEPEREL(ne)number of nodes for each element;
% ELTYPE(ne) which defines the number of element type: =0  for ragnetto
%            =3 for triangular 3 node element;
% ELSPERNODE(nn,maxelpernode)= elements around each node;
% NSTRESSPEREL(ne)= number of  stress components per element;
% NDFPEREL(ne)= number of dof per each element;
% VOL(ne) volumes of triangular plane elements;
% icrk=1 cracking modes; icrk=0 plastic mode
%
% 
            tol=10E-7
%
% Example ragnetto
%
    ragnetto=0;
     if ragnetto==1
        DIRG(1)=0
        DIRG(2)=-1
        nslc=1
        icrk=0
        ncn=5;
        ndf=1;
        prtype=1;
        nn=12;
        ne=6;
        ncoor=2;
        l=1;
        b=1;
        maxnnperel=2;
        maxelpernode=1;
        maxndfpernode=1;
        th=0
        ALFAMAX(1)=400
        gama=0
        ncrl=0
        NODESPERCRL=0;
        ONESELPERCRN=0;
        NCRNPL=0;
        RCR=0;
        TCR=0;
        INFLEN=0;
        
%
% Definition of nodal coordinates for the "ragnetto" problem
%
        COORD(1,1)=0;
        COORD(1,2)=0;
        COORD(2,1)=l;
        COORD(2,2)=0;
        COORD(3,1)=2*b;
        COORD(3,2)=0;
        COORD(4,1)=2*b;
        COORD(4,2)=l;
        COORD(5,1)=4*b;
        COORD(5,2)=0;
        COORD(6,1)=4*b;
        COORD(6,2)=l;
        COORD(7,1)=3*b;
        COORD(7,2)=l;
        COORD(8,1)=3*b;
        COORD(8,2)=2*l;
        COORD(9,1)=5*b;
        COORD(9,2)=0;
        COORD(10,1)=5*b;
        COORD(10,2)=2*l;
        COORD(11,1)=7*b;
        COORD(11,2)=0;
        COORD(12,1)=7*b;
        COORD(12,2)=2*l;
% 
% Definition of element nodes ELNODES(NE,2) 
% and element type ELTYPE(NE,1); ELTYPE=0 for the ragnetto problem
%
        ELNODES(1,1)=1;
        ELNODES(1,2)=2;
        ELTYPE(1,1)=0;
        ELNODES(2,1)=3;
        ELNODES(2,2)=4;
        ELTYPE(2,1)=0;
        ELNODES(3,1)=5;
        ELNODES(3,2)=6;
        ELTYPE(3,1)=0;
        ELNODES(4,1)=7;
        ELNODES(4,2)=8;
        ELTYPE(4,1)=0;
        ELNODES(5,1)=9;
        ELNODES(5,2)=10;
        ELTYPE(5,1)=0;
        ELNODES(6,1)=11;
        ELNODES(6,2)=12;
        ELTYPE(6,1)=0;
     end % if raagnetto==1
%
% Example 4 triangles
%
  fourtriangles=0
  if fourtriangles==1
        DIRG(1)=0
        DIRG(2)=-1
        icrk=1
        nslc=2
        ALFAMAX(1)=1
        ALFAMAX(2)=10
        ncn=7;
        ndf=2;
        prtype=6;
        nn=5;
        ne=4;
        ncor=2;
        maxnnperel=3;
        maxelpernode=4;
        th=200
       gama=2*10^-5
 %    gama=0
        ncrl=0
        NODESPERCRL=0;
        ONESELPERCRN=0;
        NCRNPL=0;
        RCR=0;
        TCR=0;
        INFLEN=0;
% Nodal coordinates
%
        COORD(1,1)=0;
        COORD(1,2)=0;
        COORD(2,1)=0;
        COORD(2,2)=200;
        COORD(3,1)=200;
        COORD(3,2)=0;
        COORD(4,1)=200;
        COORD(4,2)=200;
        COORD(5,1)=100;
        COORD(5,2)=100;
%
%  Nodes per element ELNODES(ne,3) (anticlockwise)
%
    ELNODES(1,1)=5;
    ELNODES(1,2)=2;
    ELNODES(1,3)=1;
    ELNODES(2,1)=1;
    ELNODES(2,2)=5;
    ELNODES(2,3)=3;
    ELNODES(3,1)=5;
    ELNODES(3,2)=4;
    ELNODES(3,3)=3;
    ELNODES(4,1)=2;
    ELNODES(4,2)=5;
    ELNODES(4,3)=4;
  end % if fourtriangles==1
%
% 2notch example
%
  twonotch=1
  if twonotch==1
        DIRG(1)=0
        DIRG(2)=-1
        icrk=1
        nslc=1
        ALFAMAX(1)=10
        ncn=7;
        ndf=2;
        prtype=6;
        nn=529;
        ne=956;
        ncor=2;
        maxnnperel=3;
        maxelpernode=6;
        th=50
        gama=0
        ncrl=0
        NODESPERCRL=0;
        ONESELPERCRN=0;
        NCRNPL=0;
        RCR=0;
        TCR=0;
        INFLEN=0;
%
      load ('Coordinates.txt')
        COORD=Coordinates;
%       
      load ('Element_matrix.txt')
         ELNODES=Element_matrix;
%
  end % if twonotch==1
  
          if fourtriangles ==1
            for i=1:ne
                ELTYPE(i,1)=3;
            end
          end 
%         
          if  twonotch==1
            for i=1:ne
                ELTYPE(i,1)=3;
            end
          end 
%         
%
%
% Example: colface 
%
   colface=0
   if colface==1
       load ('Element_matrix.txt');
       ELNODES=Element_matrix
       load('Coordinates.txt')
       COORD=Coordinates;
   end
%   
%
%
%  Definition of the vector NNODEPEREL(ne) which list the numeber of
%  nodes for each element; for the ragnetto problem NNODEPEREL=2; for
%  triangular 3 nodes NNODEPEREL=3;
%
        for i=1:ne
            if ELTYPE(i,1)==0
            NNODEPEREL(i,1)=2;
            end
            if ELTYPE(i,1)==3
                NNODEPEREL(i,1)=3;
            end
        end
 %
 % Definition of vector NSTRESSPEREL(ne) and tnelstress
 %
        tnelstress=0;
        for i=1:ne
            if ELTYPE(i,1)==0
                NSTRESSPEREL(i,1)=1;
                tnelstress=tnelstress+NSTRESSPEREL(i,1);
            end
            if ELTYPE(i,1)==3
                NSTRESSPEREL(i,1)=3;
                tnelstress=tnelstress+NSTRESSPEREL(i,1);
            end

        end
        
%
% computation of matrix of elements around a node ELSPERNODE(nn,maxnnperel)
%
    for i=1:nn
        r=0;
        for j=1:ne
           
            for k=1:maxnnperel;
                if ELNODES(j,k)==i;
                r=r+1;
                ELSPERNODE(i,r)=j;
                end
            end
         end
    end
%
% Definition of total number of structure dof, not accouting for
% constraints
%
        tndof=0;
        for i=1:nn;
            
            tndof=tndof+ndf;
        end
%
% Total number of elements dof tneldof
%
        tneldof=0;
        
        for i=1:ne
            nnode=NNODEPEREL(i,1);
            for j=1:nnode
                numnode=ELNODES(i,j);
                
                for k=1:ndf
                   tneldof=tneldof+1;
                end
            end
         end
    end  