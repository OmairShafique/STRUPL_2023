classdef Analysis

    properties
        number_of_dof_per_node = 0;
        dim = 0;
        Elastic_Modulus = 0;
        poissons_ratio = 0; 
        fm = 0; % compressive strenght of mortar
        ft =0; % tensile strength of mortar
        
        sigma_t = 0;
        Percentage_Limit_Tension = 0;
        tau = 0;
        Percentage_Limit_Shear = 0;
        Percentage_Fracture_Energy_Traction = 0;
        d_u = 0;
        Percentage_Ultimate_Deformation_Maximum_Step = 0;
        tol = 0;

        X_origin = 0;
        Y_origin = 0;


    end

end