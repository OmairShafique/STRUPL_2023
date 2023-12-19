function [TC,TF]=tmcs(tndof,tncdof,tnfdof,CDOF,FDOF)
%
% funcion which determines the topology matrices of the constrained 
% structure;
%
% tndof= total number of dof;
% tncdof= total number of constrained dof;
% tnfdof= total number of free dof; 
%
% TC(tndof,tncdof)= topological matrix which connect the uncontrained 
% structural dofdof to the constrained structural dof;
% TF(tndof, tnfdof)= topological matrix which connects the unconstrained
% structural dof to the free structural dof;
%
    for i=1: tndof
        for j=1:tncdof
            TC(i,j)=0;
        end
    end
    for i=1:tndof
        for j=1:tnfdof
            TF(i,j)=0;
        end
    end


    for j=1:tncdof
        indc=CDOF(j);
        TC(indc,j)=1;
    end
    for j=1:tnfdof
        indf=FDOF(j);
        TF(indf,j)=1;
    end

end
