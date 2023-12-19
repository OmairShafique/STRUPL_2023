
function [BE]=bemx(nel,nste,ndfe,elt);
if elt==0
    
    for i=1:nste
            for j=1:ndfe
                if j==1
                    BE(i,j)=-1;
                else
                    BE(i,j)=1;
                end
            end
     end
   
% if nel==5
%    for i=1:nste
%            for j=1:ndfe
%                if j==1
%                    BE(i,j)=-1;
%                else
%                    BE(i,j)=1;
%                end
%            end
%     end
% end
% if nel==6
%    for i=1:nste
%            for j=1:ndfe
%                if j==1
%                    BE(i,j)=-1;
%                else
%                    BE(i,j)=1;
%                end
%            end
%     end
%   end
end

