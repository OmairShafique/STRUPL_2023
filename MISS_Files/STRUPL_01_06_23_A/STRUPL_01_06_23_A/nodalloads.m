function [FSFM,FS]=nodalloads(COORD,ELNODES,ELTYPE,nn,ndf,ne,gama,VOL,...
    DIRG,ragnetto,TFM,TF,fourtriangles,inslc,twonotch)
%
% Computes vector FS of all (free+constrained nodal loads) 
%
%
% The case of fourtriangles plane stress problem
% 
% General case of plane 3 nodes triangular plane discretization 
%
    if ragnetto ==0
         if inslc==1
                    for i=1:nn
                        rowindex=2*i-1;
                        FS(rowindex,inslc)=0;
                        FS(rowindex+1,inslc)=0;
                    end
%
                if twonotch==0
                    [FS]=gravityloads(COORD,ELNODES,ELTYPE,nn,ndf,ne,gama,VOL,DIRG);
                end % if twonotch
%             
     if twonotch==1
         FS(2*518+2)=219.23;
         FS(2*524+2)=181.6;
                for i=1:5
                    FS(2*519+i*2-1)=0;
                    FS(2*519+i*2)=363.2;
                end %for
      end %if twonotch
          end % if inslc==1
            if inslc>1
%                
            if fourtriangles ==1
%                    
                    FS(4,inslc)=-200000;
                    FS(8,inslc)=-200000;
            end % if fourtriangles

            
%
         end %if inslc
%
    end %if ragnetto

% The case of ragnetto
%
       if ragnetto ==1
                    FS(1,inslc)=0;
                    FS(2,inslc)=0;
                    FS(3,inslc)=0;
                    FS(4,inslc)=0;
                    FS(5,inslc)=0;
                    FS(6,inslc)=0;
                    FS(7,inslc)=0;
                    FS(8,inslc)=0;
                    FS(9,inslc)=0;
                    FS(10,inslc)=0.5;
                    FS(11,inslc)=0;
                    FS(12,inslc)=0.5;
            end
        

%
%   computes nodal free loads FSF
%
            FSF=TF'*FS(:,inslc);
            
%   computes  free master nodal loads
%
            FSFM=TFM'*FSF;
%
end

    