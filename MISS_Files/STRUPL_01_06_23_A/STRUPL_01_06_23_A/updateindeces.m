function [INB,OUTB,inrow]=updateindeces(INB,OUTB,ny)
%
% display ('entering in updateindeces')
ny1=ny+1;
% inrow=index of the variable exiting out of Basis 
inrow=INB(ny);
%
% outcol=index of the variable entering into Basis
%
outcol=OUTB(ny1);
%
INB(ny)=outcol;
OUTB(ny1)=inrow;
% display ('Out from upindeces')
end
