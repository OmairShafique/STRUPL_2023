classdef Element

    %% Properties
    properties
        Length_of_Element = 0;
        Width_of_Element = 0;

        NODES;

        dhx = 0; %size in the x direction
        dhy = 0; %size in the y direction
    end

    %% Methods
    methods
        function element = Element(Length_of_Element,Width_of_Element,NODES)
            element.Length_of_Element = Length_of_Element;
            element.Width_of_Element = Width_of_Element; 
            element.NODES = NODES;
        end
    end
end

