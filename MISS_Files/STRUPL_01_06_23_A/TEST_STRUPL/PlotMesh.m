function PlotMesh(coordinates,nodes)

%%INTRODUCE 3D
%--------------------------------------------------------------------------
% Code written by : Siva Srinivas Kolukula                                |
%                   Senior Research Fellow                                |
%                   Structural Mechanics Laboratory                       |
%                   Indira Gandhi Center for Atomic Research              |
%                   India                                                 |
% E-mail : allwayzitzme@gmail.com                                         |
%          http://sites.google.com/site/kolukulasivasrinivas/             |    
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Purpose:
%         To plot the Finite Element Method Mesh
% Synopsis :
%           PlotMesh(coordinates,nodes)
% Variable Description:
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [node X Y] 
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]    
%--------------------------------------------------------------------------

nel = length(nodes) ;                  % number of elements
nnd = length(coordinates);          % total number of nodes in system
nne = size(nodes,2); 
% coordinates=geom;% number of nodes per element
%
% Initialization of the required matrices
X = zeros(nne,nel) ;
Y = zeros(nne,nel) ;

for iel=1:nel
    for i=1:nne
        nnd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
        X(i,iel)=coordinates(nnd(i),1);    % extract x value of the node
        Y(i,iel)=coordinates(nnd(i),2);    % extract y value of the node
    end
end
if nne==8
    
    % Plotting the FEM mesh, diaplay Node numbers and Element numbers
    f1 = figure ;
    set(f1,'name','Mesh','numbertitle','off') ;
    plot(X,Y,'k')
    fill(X,Y,'w')
    
    title('Finite Element Mesh') ;
    axis off ;
    k = nodes(:,1:end);
    nnd = k' ;
    for i = 1:nel
        text(X(:,i),Y(:,i),int2str(nnd(:,i)),'fontsize',8,'color','k');
        text(sum(X(:,i))/8,sum(Y(:,i))/8,int2str(i),'fontsize',10,'color','r') ;
    end
elseif nne==4
    % Plotting the FEM mesh, diaplay Node numbers and Element numbers
    f1 = figure ;
    set(f1,'name','Mesh','numbertitle','off') ;
    plot(X,Y,'k')
    fill(X,Y,'w')
    
    title('Finite Element Mesh') ;
    axis off ;
    k = nodes(:,1:end);
    nnd = k' ;
    for i = 1:nel
        text(X(:,i),Y(:,i),int2str(nnd(:,i)),'fontsize',8,'color','k');
        text(sum(X(:,i))/4,sum(Y(:,i))/4,int2str(i),'fontsize',10,'color','r') ;
    end
elseif nne==3
    % Plotting the FEM mesh, diaplay Node numbers and Element numbers
    f1 = figure ;
    set(f1,'name','Mesh','numbertitle','off') ;
    plot(X,Y,'k')
    fill(X,Y,'w')
    
    title('Finite Element Mesh') ;
    axis off ;
    k = nodes(:,1:end);
    nnd = k' ;
    for i = 1:nel
        text(X(:,i),Y(:,i),int2str(nnd(:,i)),'fontsize',8,'color','k');
        text(sum(X(:,i))/3,sum(Y(:,i))/3,int2str(i),'fontsize',10,'color','r') ;
    end
end
