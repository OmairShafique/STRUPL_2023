function [coord,g] = platelem_q4(i,analysisObject)
%
% This function returns the coordinates of the nodes of element i
% and its steering vector
% %


% To cchange the size of the problem or change the elastic properties
% ALTER the PlateQ8_input_module.m
%
%in form_KK where the KK() is trying to access a g() with a dimension = eldof but
%is restricted to only number_of_nodes_per_element or dim

coord=zeros(analysisObject.STRUCTURE.number_of_nodes_per_element,analysisObject.dim);
for k=1: analysisObject.STRUCTURE.number_of_nodes_per_element
    for j=1:analysisObject.dim
        coord(k,j)=analysisObject.geom(analysisObject.connec(i,k),j);
    end
end
%
l=0;
g = 0;
for k=1:analysisObject.STRUCTURE.number_of_nodes_per_element
    for j=1:analysisObject.number_of_dof_per_node
        l=l+1;
        g(l)=nf_g(analysisObject.connec(i,k),j);
    end
end

% End function platelem_q4