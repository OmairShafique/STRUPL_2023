function [bee,fun_3,g,A,d_3] = elem_T3(iel,analysisObject)
%
% This function returns the coordinates of the nodes of element i
% and its steering vector
%

% if Element_Type==3
    x1 = analysisObject.geom(analysisObject.connec(iel,1),1); y1 = analysisObject.geom(analysisObject.connec(iel,1),2);
    x2 = analysisObject.geom(analysisObject.connec(iel,2),1); y2 = analysisObject.geom(analysisObject.connec(iel,2),2);
    x3 = analysisObject.geom(analysisObject.connec(iel,3),1); y3 = analysisObject.geom(analysisObject.connec(iel,3),2);
    %
    %
    d_3=det([1 x1 y1; ...
        1 x2 y2; ...
        1 x3 y3]);

    A = (0.5)*det([1 x1 y1; ...
        1 x2 y2; ...
        1 x3 y3]);


    m11 = (x2*y3 - x3*y2)/(2*A);
    m21 = (x3*y1 - x1*y3)/(2*A);
    m31 = (x1*y2 - y1*x2)/(2*A);
    m12 = (y2 - y3)/(2*A);
    m22 = (y3 - y1)/(2*A);
    m32 = (y1 - y2)/(2*A);
    m13 = (x3 - x2)/(2*A);
    m23 = (x1 - x3)/(2*A);
    m33 = (x2 -x1)/(2*A);

    bee = [ m12 0 m22 0 m32 0; ...
        0 m13 0 m23 0 m33; ...
        m13 m12 m23 m22 m33 m32] ;

    fun_3= [(m11+m12+m13)  0;...
        0  (m11+m12+m13);...
        (m21+m22+m23)  0;...
        0  (m21+m22+m23);...
        (m31+m32+m33)  0;...
        0  (m31+m32+m33)];
    %

    l=0;
    for k=1: analysisObject.STRUCTURE.number_of_nodes_per_element
        for j=1:analysisObject.number_of_dof_per_node
            l=l+1;
            g(l)=analysisObject.nf_g(analysisObject.connec(iel,k),j);
        end
    end
