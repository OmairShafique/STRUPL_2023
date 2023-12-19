function [TA,INB,OUTB]=rearrange(TA,INB,OUTB,indexi,indexj,ny);
%
% Store temporary row i,row ny, column j and column 2+ny
% into four temporary vectors; then exchange row i with row ny and column 
% j with column 2+ny in matrix TA and update vectors INB and OUTB
%
%
% Rearrange 2 rows
%
% display ('entro in rearrange')
%
rowi=TA(indexi,:);
rowny=TA(ny,:);
TA(indexi,:)=rowny;
TA(ny,:)=rowi;
%
% Rearrange 2 columns
%
columnj=TA(:,indexj);
ny2=ny+2;
columnny2=TA(:,ny2);
TA(:,indexj)=columnny2;
TA(:,ny2)=columnj;
%
%
% rearrange vectors INB and OUTB
% 
% vector INB
%
%
iout=INB(indexi,1);
iny=INB(ny);
INB(indexi,1)=iny;
INB(ny,1)=iout;
%
% vector OUTB
%
indexj1=indexj-1;
ny1=ny+1;
jin=OUTB(1,indexj1);
jny=OUTB(1,ny1);
OUTB(1,ny1)=jin;
OUTB(1,indexj1)=jny;
%
% display ('esco da rearrange')
end
