function [indexj]=inbasis(OUTB,ny,inrow);
% if inrow <=NY then incol=inrow+NY and INDEXJ=k where k is such that
% OUTB(k)=incol
% if inrow >NY then incol= 2*NY-inrow and INDEXJ=k as before.
%
% display ('entro in inbasis')
if inrow<= ny
    incol=inrow+ny;
else
    incol=inrow-ny;
end
ny1=ny+1;
for k=1:ny1
    if OUTB(k)==incol;
indexj=k+1;
    end
% display ('esco da inbasis')
end
