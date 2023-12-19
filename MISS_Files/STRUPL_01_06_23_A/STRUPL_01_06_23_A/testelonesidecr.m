%
% Program for testing the function "elspercrnode" which defines the element
% numbers along one side of a linear crack;
%
% definition of the nodal coordinates of triangular elemnts around one node
%
% coordinates of node i=1
    COORD(1,1)=1;
    COORD(1,2)=1;
%
% coordinates of nodes around node i; nodes 2,3,4,5,6,7,8,9 matrix COORD
%
    COORD(2,1)=0;
    COORD(2,2)=0;
    COORD(3,1)=1;
    COORD(3,2)=0;
    COORD(4,1)=2;
    COORD(4,2)=0;
    COORD(5,1)=2;
    COORD(5,2)=1;
    COORD(6,1)=2;
    COORD(6,2)=2;
    COORD(7,1)=1;
    COORD(7,2)=2;
    COORD(8,1)=0;
    COORD(8,2)=2;
    COORD(9,1)=0;
    COORD(9,2)=1;
%
% definition of the nodes per each element:matrix ELSNODES
%
    ELNODES(1,1)=1;
    ELNODES(1,2)=5;
    ELNODES(1,3)=6;
    ELNODES(2,1)=1;
    ELNODES(2,2)=6;
    ELNODES(2,3)=7;
    ELNODES(3,1)=1;
    ELNODES(3,2)=7;
    ELNODES(3,3)=8;
    ELNODES(4,1)=9;
    ELNODES(4,2)=1;
    ELNODES(4,3)=8;
    ELNODES(5,1)=2;
    ELNODES(5,2)=1;
    ELNODES(5,3)=9;
    ELNODES(6,1)=2;
    ELNODES(6,2)=3;
    ELNODES(6,3)=1;
    ELNODES(7,1)=3;
    ELNODES(7,2)=4;
    ELNODES(7,3)=1;
    ELNODES(8,1)=4;
    ELNODES(8,2)=5;
    ELNODES(8,3)=1;
    
%
% Definition of elements aroud a node 1 ELSPERNODE
%
    ELSPERNODE(1,1)= 1;
    ELSPERNODE(1,2)= 2;
    ELSPERNODE(1,3)= 3;
    ELSPERNODE(1,4)= 4;
    ELSPERNODE(1,5)= 5;
    ELSPERNODE(1,6)= 6;
    ELSPERNODE(1,7)= 7;
    ELSPERNODE(1,8)= 8;
    ELSPERNODE(2,1)= 5;
    ELSPERNODE(2,2)= 6;
    ELSPERNODE(2,3)= 0;
    ELSPERNODE(2,4)= 0;
    ELSPERNODE(2,5)= 0;
    ELSPERNODE(2,6)= 0;
    ELSPERNODE(2,7)= 0;
    ELSPERNODE(2,8)= 0;
    ELSPERNODE(3,1)= 6;
    ELSPERNODE(3,2)= 7;
    ELSPERNODE(3,3)= 0;
    ELSPERNODE(3,4)= 0;
    ELSPERNODE(3,5)= 0;
    ELSPERNODE(3,6)= 0;
    ELSPERNODE(3,7)= 0;
    ELSPERNODE(3,8)= 0;
    ELSPERNODE(4,1)= 7;
    ELSPERNODE(4,2)= 8;
    ELSPERNODE(4,3)= 0;
    ELSPERNODE(4,4)= 0;
    ELSPERNODE(4,5)= 0;
    ELSPERNODE(4,6)= 0;
    ELSPERNODE(4,7)= 0;
    ELSPERNODE(4,8)= 0;
    ELSPERNODE(5,1)= 8;
    ELSPERNODE(5,2)= 1;
    ELSPERNODE(5,3)= 0;
    ELSPERNODE(5,4)= 0;
    ELSPERNODE(5,5)= 0;
    ELSPERNODE(5,6)= 0;
    ELSPERNODE(5,7)= 0;
    ELSPERNODE(5,8)= 0;
    ELSPERNODE(6,1)= 1;
    ELSPERNODE(6,2)= 2;
    ELSPERNODE(6,3)= 0;
    ELSPERNODE(6,4)= 0;
    ELSPERNODE(6,5)= 0;
    ELSPERNODE(6,6)= 0;
    ELSPERNODE(6,7)= 0;
    ELSPERNODE(6,8)= 0;
    ELSPERNODE(7,1)= 2;
    ELSPERNODE(7,2)= 3;
    ELSPERNODE(7,3)= 0;
    ELSPERNODE(7,4)= 0;
    ELSPERNODE(7,5)= 0;
    ELSPERNODE(7,6)= 0;
    ELSPERNODE(7,7)= 0;
    ELSPERNODE(7,8)= 0;
    ELSPERNODE(8,1)= 3;
    ELSPERNODE(8,2)= 4;
    ELSPERNODE(8,3)= 0;
    ELSPERNODE(8,4)= 0;
    ELSPERNODE(8,5)= 0;
    ELSPERNODE(8,6)= 0;
    ELSPERNODE(8,7)= 0;
    ELSPERNODE(8,8)= 0;
    ELSPERNODE(9,1)= 4;
    ELSPERNODE(9,2)= 5;
    ELSPERNODE(9,3)= 0;
    ELSPERNODE(9,4)= 0;
    ELSPERNODE(9,5)= 0;
    ELSPERNODE(9,6)= 0;
    ELSPERNODE(9,7)= 0;
    ELSPERNODE(9,8)= 0;
%
     maxelspernode=8
%
% direction of the normal to the crack directions
%
      ncrl=4
      BETA(3,1)=  0*pi/180
      BETA(2,1)= 45*pi/180
      BETA(1,1)= -90*pi/180
      BETA(4,1)=  -45*pi/180
%      
%    Vector NCRNPL; vector which gives the number of nodes per each
%    cracking line
%
        NCRNPL(1)=3
        NCRNPL(2)=3
        NCRNPL(3)=3
        NCRNPL(4)=3
       
% Matrix NODESPERCRL
%
        NODESPERCRL(1,1)=2
        NODESPERCRL(1,2)=3
        NODESPERCRL(1,3)=4
        NODESPERCRL(2,1)=4
        NODESPERCRL(2,2)=1
        NODESPERCRL(2,3)=8
        NODESPERCRL(3,1)=3
        NODESPERCRL(3,2)=1
        NODESPERCRL(3,3)=7
        NODESPERCRL(4,1)=2
        NODESPERCRL(4,2)=1
        NODESPERCRL(4,3)=6
%
%
% Function elpercrnode definitioni of the elements around node i
% on one side of the crack with inclinatiion beta of the normal to the 
% crack line with respect axis x1;
%
      [ONESELPERCRN]=oneselpercrn(ncrl,NCRNPL,NODESPERCRL,COORD,...
              ELNODES,ELSPERNODE,BETA,maxelspernode)
%
%
    return