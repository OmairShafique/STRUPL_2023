function [K,f] = assem(nne,ndof,C,K,Ke,f,fe)
% ASSEM -->  Assemble element matrices Ke ( and fe ) into the global
%  stiffness matrix K ( and the global force vector f )
%  according to the topology matrix edof.
%-------------------------------------------------------------
% INPUT:    nne  : Number of element nodes
%           ndof : Number of node dof 
%           C    : Element matrix
%           K    : Global stiffness matrix
%           Ke   : Element stiffness matrix
%           f    : Global force vector
%           fe   : Element force vector
%
% OUTPUT:   K    : New global stiffness matrix
%           f    : New global force vector
% -------------------------------------------------------------

% INPUT DATA -------------------------------------------------

    for aa = 1 : nne %Number of element nodes
        for ii = 1 : ndof %Number of node dof 
            rw = ndof*(C(aa)-1)+ii ; %Row position
            f(rw) = f(rw) + fe(ndof*(aa-1)+ii) ;
            for bb = 1 : nne %Number of element nodes
                for kk = 1 : ndof %Number of node dof  
                     cl = ndof*(C(bb)-1)+kk ; %Column position
                     K(rw,cl) = K(rw,cl) + Ke(ndof*(aa-1)+ii,ndof*(bb-1)+kk);
                end
            end
        end
    end        
end
%--------------------------end--------------------------------

