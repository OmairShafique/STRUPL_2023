function[Global_stiffness_matrix]=form_KK(Global_stiffness_matrix, kg, g)
%
% This function assembles the global stiffness matrix
%
global Degrees_of_Freedom_Per_Element
% Degrees_of_Freedom_Per_Elementdof = 12;
%
% This function assembles the global stiffness matrix
%

for i=1:Degrees_of_Freedom_Per_Element
    if g(i) ~= 0
        for j=1: Degrees_of_Freedom_Per_Element
            if g(j) ~= 0
                Global_stiffness_matrix(g(i),g(j))= Global_stiffness_matrix(g(i),g(j)) + kg(i,j);
            end
        end
    end
end
end
%
%%%%%%%%%%%%% end function form_KK %%%%%%%%%%%%%%%%%