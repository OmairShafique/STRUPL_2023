%
% Script to test function topcrack
%
%
% Input:
    ne=4
    ndf=2;
    ncrl=3;
    NCRPL(1)=2;
    NCRPL(2)=3;
    NCRPL(3)=3;
    NODESPERCRL(1,1)=1;
    NODESPERCRL(1,2)=3;
    NODESPERCRL(2,1)=1;
    NODESPERCRL(2,2)=5;
    NODESPERCRL(2,3)=4
    NODESPERCRL(3,1)=3;
    NODESPERCRL(3,2)=5;
    NODESPERCRL(3,3)=2;
 %   
    ELNODES(1,1)=1;
    ELNODES(1,2)=2;
    ELNODES(1,3)=5;
    ELNODES(2,1)=1;
    ELNODES(2,2)=3;
    ELNODES(2,3)=5;
    ELNODES(3,1)=5;
    ELNODES(3,2)=3;
    ELNODES(3,3)=4;
    ELNODES(4,1)=5;
    ELNODES(4,2)=4;
    ELNODES(4,3)=2;
%    
    [TCR] = topcrack(ne,ndf,ncrl,NCRPL,NODESPERCRL,ELNODES)
%
return