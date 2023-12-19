function [B] = assem_B(nne,ndof,C,B,Be)
% ASSEM -->  Assemble element matrices Ke ( and fe ) into the global
%  stiffness matrix K ( and the global force vector f )
%  according to the topology matrix edof.
%-------------------------------------------------------------
% INPUT:    nne  : Number of element nodes
%           ndof : Number of node dof 
%           C    : Element matrix
%           B    : Global compatibility matrix
%           be   : Element compatibility matrix
%
% OUTPUT:   B    : New compatibility matrix
% -------------------------------------------------------------

% INPUT DATA -------------------------------------------------

    for ii = 1 : 3 %Number of deformation components 
        rw = ii ; %Row position
        for bb = 1 : nne %Number of element nodes
            for kk = 1 : ndof %Number of node dof  
                cl = ndof*(C(bb)-1)+kk ; %Column position
                B(rw,cl) = B(rw,cl) + Be(rw,ndof*(bb-1)+kk) ;
            end
        end
    end        
end
%--------------------------end--------------------------------

