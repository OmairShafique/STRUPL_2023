% Program for testing 2 functions: consdof and tvs
% function consdof creats  vector CDOF(tncdof) which list, in progressive
% the number of constrained dof.
%
nn=4
ndf=2
ncn=2
for i=1:2
    for j=1:15
        SC(i,j)=0;
    end
end
SC(1,1)=1;
SC(1,2)=1;
SC(1,3)=1;
SC(2,1)=3;
SC(2,3)=1;
%
[CDOF,FDOF,tndof,tnfdof,tncdof]=consdof(nn,ndf,ncn,SC)
%
[TC,TF]=tmcs(tndof,tncdof,tnfdof,CDOF,FDOF)
%
return