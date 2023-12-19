function [FS,alfamax]=loads(tnfdof,ragnetto)
%
% Function loads computes the nodal force vector FS(tnsdof)
%
if ragnetto==1
    alfamax=400;
        for i=1:tnfdof
            FS(i,1)=0.;
            if i==3
                FS(3,1)=1;
            end
        end
end
end

