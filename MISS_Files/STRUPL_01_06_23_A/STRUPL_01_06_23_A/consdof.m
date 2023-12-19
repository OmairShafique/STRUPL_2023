function [CDOF,FDOF,tndof,tnfdof,tncdof,UC]=consdof(nn,ndf,ncn,SC)
%
% function consdof determines vectors CDOF (constraint dof) which list the
% constraint dof indeces and FDOF (free dof) which lists the free dof
% indeces;
%
% nn=number of nodes (input);
% ndf= number od nodal dof (input);
% ncn= number of constrained nodes (input);
% SC(ncn,15)= matrix of constraints (input);
%   SC(i,1)=node i number
%   SC(i,1-2....7)=6 nodal dof;1=constrained dof =0 free dof;
%   SC(i,8-9)=two angles of the local ref system + anticlockwise;
%   SC(i,10..-15)=prescribed nodal displacements;
%
% icdof= progressive index of constrained dof;
% ifdof= progressive index of free dof;
% idof=progressive index of total dof;
%
icdof=0;
ifdof=0;
idof=0;
for i=1:nn
  for j=1:ndf
    flag=0;
    idof=idof+1;
      for k=1:ncn
        if flag==0
           icn=SC(k,1);
            if icn==i
                
                    if flag==0
                        if SC(k,1+j)==1
                            flag=1;
                            icdof=icdof+1;
                            CDOF(icdof)=idof;
                            UC(icdof,1)=SC(k,9+j);
                         end
                    end
                
             end
         end
     end
   if flag==0
      ifdof=ifdof+1;
      FDOF(ifdof)=idof;
    end
  end
end
 tncdof=icdof;
 tndof=idof;
 tnfdof=tndof-tncdof;
 end
        