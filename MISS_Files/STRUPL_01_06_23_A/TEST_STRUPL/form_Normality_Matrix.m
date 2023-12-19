function [Normality_Matrix,Hard_Soft_Matrix]=form_Normality_Matrix(crack_nodes_lines,Number_of_crack_lines)

Normality_Matrix=[];
Hard_Soft_Matrix=[];

%
irow=0;
jcol=0;
inot=0;
nno=0;
NOT=0;
coulomb=1;
notension=0;
if twonotch==1
    coulomb=0;
    notension=1
end
if coulomb==1
    %
    % Coulomb friction angle as radians
    %
    % definition of matrix N'(ny,tncrn*ndf) for Coulomb + tension cut off
    % according to normality rule in the local reference system + hardening
    % matrix;
    %
    theta =0.523598;
    for i=1:Number_of_crack_lines
        for ii=1: crack_nodes_lines(i)
            irow=irow+1;
            jcol=jcol+1;
            %     inot=inot+1;
            % mode 1 cut off
            Normality_Matrix(irow,jcol)=1;
            Normality_Matrix(irow,jcol+1)=0;
            % mode 1 twin
            %     NT(irow+1,jcol)=1;
            %     NT(irow+1,jcol+1)=0;
            %
            %       NOT(inot,1)=irow;
            %       NOT(inot,2)=irow+1;
            %
            % mode 2
            %
            Normality_Matrix(irow+1,jcol)=sin(theta);
            Normality_Matrix(irow+1,jcol+1)=cos(theta);
            %     NT(irow+2,jcol)=sin(theta);
            %     NT(irow+2,jcol+1)=cos(theta);
            % mode 2 twin
            %     NT(irow+3,jcol)=sin(theta);
            %     NT(irow+3,jcol+1)=cos(theta);
            %
            %     NOT(inot+1,1)=irow+2;
            %     NOT(inot+1,2)=irow+3;
            % mode 3
            Normality_Matrix(irow+2,jcol)=sin(theta);
            Normality_Matrix(irow+2,jcol+1)=-cos(theta);
            %      NT(irow+4,jcol)=sin(theta);
            %      NT(irow+4,jcol+1)=-cos(theta);
            % mode 3 twin
            %     NT(irow+5,jcol)=sin(theta);
            %     NT(irow+5,jcol+1)=-cos(theta);
            %
            %     NOT(inot+2,1)=irow+4;
            %     NOT(inot+2,2)=irow+5;
            %
            % mode 1
            Hard_Soft_Matrix(irow,irow)=0;
            %    HS(irow,irow)=-1000;
            % mode 1 twin
            %     HS(irow+1,irow+1)=0;
            % mode 2
            Hard_Soft_Matrix(irow+1,irow+1)=0;
            %      HS(irow+2,irow+2)=-1000;
            % mode 2 twin
            %     HS(irow+3,irow+3)=0;
            % mode 3
            Hard_Soft_Matrix(irow+2,irow+2)=0;
            %      HS(irow+4,irow+4)=-1000;
            % mode 3 twin
            %      HS(irow+5,irow+5)=0;
            %
            %      inot=inot+2;
            jcol=jcol+1;
            %      irow=irow+5;
            irow=irow+2;
        end %for ii
    end  % for i
    ny=irow;
    %   nno=ny/2;
end %if coulomb
%
% End Coulomb material and start no tension only
%
if notension==1
    irow=0;
    jcol=0;
    inot=0;
    for i=1:ncrl
        for ii=1: crack_nodes_lines(i)
            irow=irow+1;
            jcol=jcol+1;
            inot=inot+1;
            %
            % mode 1: onnly tension
            %
            Normality_Matrix(irow,jcol)=1;
            Normality_Matrix(irow,jcol+1)=0;
            % mode 1 twin
            Normality_Matrix(irow+1,jcol)=1;
            Normality_Matrix(irow+1,jcol+1)=0;
            %
            NOT(inot,1)=irow;
            NOT(inot,2)=irow+1;
            %
            %
            % mode 1
            %     HS(irow,irow)=0;
            Hard_Soft_Matrix(irow,irow)=-40476;
            if irow==1
                Hard_Soft_Matrix(1,1)=-20238;
            end

            if irow==41
                Hard_Soft_Matrix(41,41)=-20238;
            end
            % mode 1 twin
            Hard_Soft_Matrix(irow+1,irow+1)=0;
            %
            jcol=jcol+1;
            irow=irow+1;
        end %for ii
    end  % for i
    ny=irow;
    nno=ny/2;
end %if notension

%
% end no tension material
%
end %if icrack
for j=1:ny
    FI(j)=0;
end
%
end
