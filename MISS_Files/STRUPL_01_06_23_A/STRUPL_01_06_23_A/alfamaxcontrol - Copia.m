function [ALFA]=alfamaxcontrol(NY,INB,TA,ALFAMAX)
%
display ('In alfamaxcontrol')
%Control on maximum load factor ALFA<ALFAMAX
%
NY2=NY+2
i=1
while INB(i)>0
i=i+1
ALFA=TA(i,1)
INDEXALFA=i
end 
%
if ALFA>ALFAMAX
    BETA=(ALFAMAX-ALFA)/TA(INDEXALFA,NY+2)
%
% Computation of the basic vector BSOLMAX for ALFA=ALFAMAX
%
display ('Basic solution for ALFA=ALFAMAX')
display ('Active modes')
ACTIVE_INDECES=INB
BSOLMAX= TA(1:NY,1)+BETA*TA(1:NY,NY2)
return
%
display ('Out alfamaxcontrol')
end