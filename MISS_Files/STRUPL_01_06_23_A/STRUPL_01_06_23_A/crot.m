%
% script to test function crlines
%
    fourtriangles=1;
    COORD(1,1)=0;
    COORD(1,2)=0;
    COORD(2,1)=0;
    COORD(2,2)=2;
    COORD(3,1)=2;
    COORD(3,2)=0;
    COORD(4,1)=2;
    COORD(4,2)=2;
    COORD(5,1)=1;
    COORD(5,2)=1;
[ncrl,tncrn,NCRNPL,NODESPERCRL,BETA,RCRS]=crlines(fourtriangles,COORD)
return