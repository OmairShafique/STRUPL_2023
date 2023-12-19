function[MX, MY, MXY, QX, QY]=Forces_at_nodes_plate(Element_Forces)
%
% This function averages the stresses at the nodes
%
global Number_of_Nodes Number_of_Elements number_of_nodes_per_element connec
%
for k = 1:Number_of_Nodes
mx = 0. ; my = 0.; mxy = 0.; qx = 0.; qy = 0.;
ne = 0;
for iel = 1:Number_of_Elements;
for jel=1:number_of_nodes_per_element;
if connec(iel,jel) == k;
ne=ne+1;
mx = mx + Element_Forces(iel,1);
my = my + Element_Forces(iel,2);
mxy = mxy + Element_Forces(iel,3);
qx = qx + Element_Forces(iel,4);
qy = qy + Element_Forces(iel,5);
end
end
end
MX(k,1) = mx/ne;
MY(k,1) = my/ne;
MXY(k,1) = mxy/ne;
QX(k,1) = qx/ne;
QY(k,1) = qy/ne;
end