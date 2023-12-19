function[beeb] = formbeeb(deriv,number_of_nodes_per_element,Degrees_of_Freedom_Per_Element)
%
global number_of_dof_per_node
% This function assembles the matrix [beeb] from the
% derivatives of the shape functions in global coordinates
% for a thick plate element (bending action)
%
beeb=zeros(3,Degrees_of_Freedom_Per_Element);
for m=1:number_of_nodes_per_element
    k=number_of_dof_per_node*m;
    j=k-1;
    x=deriv(1,m);
    beeb(1,j)=x;
    beeb(3,k)=x;
    y=deriv(2,m);
    beeb(2,k)=y;
    beeb(3,j)=y;
end
%
% End function formbeeb