function[VS,tnconsdof]=constraints(nn,ndf,tnuncsdof,ragnetto,...
    twotriangles)
%
% display ('enter in constraints')
% It generates matrix VMX(tnuncdof, tncondof) where
% tnuncsdof=total number of uncontrained structuredof
% tnconsdof= total number of constrained dof
%
% computationo of tnuncsdof 
%
%      
% Ragnetto problem:
% Total number of unconstrained structure dof=12
% Total number of constraints =8
% Total numberof constraint structure dof=8
% Total number of structure dof=4
%
    if ragnetto==1
        tnconsdof=8
        
    end
    if twotriangles==1
        tnconnodes=2
        tncons=3
%
%
% Vector  CONDOF(tnconstraints) gives the number of constraint dof 
%
           CONDOF(1)=1
           CONDOF(2)=2
           CONDOF(3)=6
 %
 % matrix NODEDOF list in first column the nodes with at least one free dof
 % in the following columns 0 means free dof and 1 constrained;

       
     

        tnconsdof=tnuncsdof-tncons
%
% Initialization matrix VS
%
        for i=1:tnuncsdof
            for j=1:tnuncsdof;
                PVS(i,j)=0;
                PVS(i,i)=1;
            end
        end
    end

%
%  Computes matrix VS
%
    
    
    
    
    
                  
           
                    
%            
%        
%    
%
%
    if ragnetto==1
%    
%Ragnetto example
% dof 1,3,5,9,11 are fixed therefore the corresponding lines are alla zero
% i.e. do not require any intervention;
% dof 2,4,6 are slave of master dof 1,2: t e coefficient woth respwct
% column 1 is 1 while with respect 2 is given by the distance of the slave
% dof with the dof 2, assuming positive 2 if anticlockwise;
%a nalogous reasonong is repeated for slave dof 8,9,10 with respect master 3
% and 4.
%
    VS(2,1)=1;
    VS(4,1)=1;
    VS(6,1)=1;
    VS(2,2)=3;
    VS(4,2)=1;
    VS(6,2)=-1;
    VS(7,1)=1;
    VS(8,3)=1;
    VS(10,3)=1;
    VS(12,3)=1;
    VS(8,4)=3;
    VS(10,4)=1;
    VS(12,4)=-1;
   end
%
%
%       
%
% display ('exit from constraints')
end











