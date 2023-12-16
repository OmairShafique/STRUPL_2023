classdef Element

    %% Properties
    properties
        Length_of_Element = 0;
        Width_of_Element = 0;
        number_of_degree_of_freedom_per_element = 0;
        total_number_of_elements_connected = 0;
        number_of_nodes_in_the_element = 0;

        NODES;

        dhx = 0; %size in the x direction
        dhy = 0; %size in the y direction
    end

    %% Methods
    methods
        function element = Element(Length_of_Element,Width_of_Element,number_of_degree_of_freedom_per_element,total_number_of_elements_connected,number_of_nodes_in_the_element)
            element.Length_of_Element = Length_of_Element;
            element.Width_of_Element = Width_of_Element;
            element.number_of_degree_of_freedom_per_element = number_of_degree_of_freedom_per_element;
            element.total_number_of_elements_connected = total_number_of_elements_connected;
            element.number_of_nodes_in_the_element = number_of_nodes_in_the_element;
        end
    end
end

