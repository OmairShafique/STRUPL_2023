
function[ndf,COORD,ELNODES,prtype,nn,ne,maxnnperel,maxndfpernode, ...
    maxelpernode,ELSPERNODE,NNODEPEREL,tneldof,tnuncsdof,ELTYPE,...
    NSTRESSPEREL,NDFPEREL,ragnetto,twotriangles]=geometry;
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
%
% ndf= number of nodal degrees of freedom n.d.f.
%   ndf=2 for plane truss and plane stress and strain problems ;
%   ndf=3 for plane frames, space brick and plate problems ;
%   ndf=6 for space frame and shell problems
%
% maxnnperel= max number of nodes per elements
% tneldof = Total number of element dof
% tnuncsdof= total number of unconstrained structure dof
% tnelstress= total number of elemnt stress components
% maxndfpernode= Max number of dof per node
% maxelpernode= Max number of elements around a node
% maxnnperel= Max number of nodes per element
%
%             Vectors and matrices
% COORD(nn,ncoor)  nodal coordinates;
% ELNODES(ne, maxnnperel)nodes of each  element;
% ELSPERNODE (nn,maxelpernode) elements around each node;
% NNODEPEREL(ne)number of nodes for each elemnt;
% 
% ELTYPE(ne) which defines the number of element type: =0  for ragnetto
%            =3 for triangular 3 node element;
% ELSPERNODE(nn,maxelpernode)= elements around each node;
% NSTRESSPEREL(ne)= number of  stress components per element;
% NDFPEREL(ne)= number of dof per each element;
%
%
% 
    ragnetto=0;
   if ragnetto == 1
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
% and element type ELTYPE(NE); ELTYPE=0 for the ragnetto problem
%
        ELNODES(1,1)=1;
        ELNODES(1,2)=2;
        ELTYPE(1)=0;
        ELNODES(2,1)=3;
        ELNODES(2,2)=4;
        ELTYPE(2)=0;
        ELNODES(3,1)=5;
        ELNODES(3,2)=6;
        ELTYPE(3)=0;
        ELNODES(4,1)=7;
        ELNODES(4,2)=8;
        ELTYPE(4)=0;
        ELNODES(5,1)=9;
        ELNODES(5,2)=10;
        ELTYPE(5)=0;
        ELNODES(6,1)=11;
        ELNODES(6,2)=12;
        ELTYPE(6)=0;
   end
%
%   Data of problem 2triangles
%
    twotriangles=1
    if twotriangles==1;
        ndf=2
        prtype=2
        nn=4
        ne=2
        ncor=2
        maxnnperel=3
        maxelpernode=2
        maxndfpernode=2
%
% Nodal coordinates
%
        COORD(1,1)=0
        
        COORD(1,2)=0
        COORD(2,1)=4
        COORD(2,2)=0
        COORD(3,1)=2
        COORD(3,2)=0
        COORD(4,1)=2
        COORD(4,2)=4
%
%  Elements per node ELNODES(ne)
%
    ELNODES(1,1)=1;
    ELNODES(1,2)=2;
    ELNODES(1,3)=4;
    ELNODES(2,1)=1;
    ELNODES(2,2)=3;
    ELNODES(2,3)=4;
%
%
%   Element type
%
       for i=1:ne
           ELTYPE(i)=3;
       end
%
    end

%
%  Definition of the vector NNODESPERELEMENT(ne) which list the numeber of
%  nodes for each element; for the ragnetto problem NNODEPEREL=2; for
%  triangular 3 nodes NNODEPEREL=3;
%
        for i=1:ne
            if ELTYPE(i)==0
            NNODEPEREL(i)=2;
            end
            if ELTYPE(i)==3
                NNODEPEREL(i)=3;
            end
        end
 %
 % Definition of vector NSTRESSPEREL(ne) and tnelstress
 %
        tnelstress=0;
        for i=1:ne
            if ELTYPE(i)==0
                NSTRESSPEREL(i)=1;
                tnelstress=tnelstress+NSTRESSPEREL(i);
            end
            if ELTYPE(i)==3
                NSTRESSPEREL(i)=3
                tnelstress=tnelstress+NSTRESSPEREL(i)
            end

        end
        
%
% computation of matrix of elements around a node ELSPERNODE(nn,maxnnperel)
%
    for i=1:nn
        r=0
        for j=1:ne
           
            for k=1:maxnnperel;
                if ELNODES(j,k)==i;
                r=r+1
                ELSPERNODE(i,r)=j
                end
            end
         end
    end
%
% Definition of total number of structure dof, not accouting for
% constraints
%
        tnuncsdof=0;
        for i=1:nn;
            
            tnuncsdof=tnuncsdof+ndf
        end
        
        
%
% Total number of elements dof tneldof
%
        tneldof=0
        for i=1:ne
            NDFPEREL(i)=0
        end
        for i=1:ne
            nnode=NNODEPEREL(i)
            for j=1:nnode
                numnode=ELNODES(i,j)
                
                NDFPEREL(i)=NDFPEREL(i)+ndf
                for k=1:ndf
                   tneldof=tneldof+1
                end
            end
         end
  
  
%
return
    end