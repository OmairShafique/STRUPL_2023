function[coord,g] = platelem_q4(i)
%
% This function returns the coordinates of the nodes of element i
% and its steering vector
% %

global number_of_nodes_per_element number_of_dof_per_node geom connec nf_g  dim
% To cchange the size of the problem or change the elastic properties
% ALTER the PlateQ8_input_module.m
%
%in form_KK where the KK() is trying to access a g() with a dimension = eldof but
%is restricted to only number_of_nodes_per_element or dim

coord=zeros(number_of_nodes_per_element,dim);
for k=1: number_of_nodes_per_element
    for j=1:dim
        coord(k,j)=geom(connec(i,k),j);
    end
end
%
l=0;
g = 0;
for k=1:number_of_nodes_per_element
    for j=1:number_of_dof_per_node
        l=l+1;
        g(l)=nf_g(connec(i,k),j);
    end
end

% End function platelem_q4