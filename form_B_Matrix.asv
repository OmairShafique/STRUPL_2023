function[B_Matrix] = form_B_Matrix(number_of_nodes_per_element,Degrees_of_Freedom_Per_Element,Element_type,ngpb,ngps)
% This function assembles the matrix B for plane problems from the
% derivatives of the shape functions in global coordinates
% for a thick plate element (bending action)
global number_of_nodes_per_element Degrees_of_Freedom_Per_Element Element_type ngpb ngps geom connec


if Element_type==3 && ngpb==0

    x1 = geom(connec(iel,1),1); y1 = geom(connec(iel,1),2);
    x2 = geom(connec(iel,2),1); y2 = geom(connec(iel,2),2);
    x3 = geom(connec(iel,3),1); y3 = geom(connec(iel,3),2);
    %
    A = (0.5)*det([1 x1 y1; ...
        1 x2 y2; ...
        1 x3 y3]);
    %
    m11 = (x2*y3 - x3*y2)/(2*A);
    m21 = (x3*y1 - x1*y3)/(2*A);
    m31 = (x1*y2 - y1*x2)/(2*A);
    m12 = (y2 - y3)/(2*A);
    m22 = (y3 - y1)/(2*A);
    m32 = (y1 - y2)/(2*A);
    m13 = (x3 - x2)/(2*A);
    m23 = (x1 - x3)/(2*A);
    m33 = (x2 -x1)/(2*A);
    %
    B_Matrix = [ m12 0 m22 0 m32 0; ...
        0 m13 0 m23 0 m33; ...
        m13 m12 m23 m22 m33 m32] ;
end
if Element_type~=3 && ngpb~=0

    B_Matrix=zeros(3,Degrees_of_Freedom_Per_Element);
    for m=1:number_of_nodes_per_element
        k=3*m;
        j=k-1;
        x=deriv(1,m);
        B_Matrix(1,j)=x;
        B_Matrix(3,k)=x;
        y=deriv(2,m);
        B_Matrix(2,k)=y;
        B_Matrix(3,j)=y;
    end
end
if Element_type~=3 && ngps~=0
    B_Matrix=zeros(2,Degrees_of_Freedom_Per_Element);
    for m=1:number_of_nodes_per_element
        k=3*m;
        j=k-1;
        i=k-2;
        x=deriv(1,m); y=deriv(2,m);
        B_Matrix(2,i)=-x;
        B_Matrix(1,i)=-y;
        B_Matrix(1,k) = fun(m);
        B_Matrix(2,j) = fun(m);
    end
end
% End function form_B_Matrix