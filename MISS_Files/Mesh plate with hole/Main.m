% To mesh a plate with hole at center using Transfinite Interpolation (TFI)
%{
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Warning : On running this the workspace memory will be deleted. Save if
 any data present before running the code !!
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
--------------------------------------------------------------------------
 Code written by : Siva Srinivas Kolukula, PhD                           |
                   Structural Mechanics Laboratory                       |
                   Indira Gandhi Center for Atomic Research              |
                   India                                                 |
 E-mail : allwayzitzme@gmail.com                                         |
 web-link: https://sites.google.com/site/kolukulasivasrinivas/           |
--------------------------------------------------------------------------
%}
% Version 1 : 15 November 2013
clc ; clear all
% Dimensions of the plate
L = 1. ;                % Length of the plate
B = 1. ;                % Breadth of the plate

% Number of discretizations along xi and eta axis
m = 10 ;
n = 10 ;
%
% Model plate as two regions which lie in first quadrant
global R theta;
R = 0.1 ;               % Radius of the hole at center
%%%%%%%%%%%%%%%%%%%%%%%Dont change from here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = pi/2 ;          % Quarter angle of the hole
global O P1 P2 P3 P4 P5 CMP ;
O = [0. 0.] ;           % Centre of plate and hole
P1 = [R 0.] ;           % Edge of the hole and plate
P2 = [L/2 0.] ;         % Edge of the plate
P3 = [L/2 B/2] ;        % Edge of the plate
P4 = [0. B/2] ;         % Edge of the plate
P5 = [0. R] ;           % Edge of the hole and plate
CMP = [R*cos(theta/2.) R*sin(theta/2.)] ;
% discretize along xi and eta axis
xi = linspace(0.,1,m) ;
eta = linspace(0.,1.,n) ;
% Number of Domains 
Domain = 2 ;
DX = cell(1,Domain) ;   
DY = cell(1,Domain) ;
for d = 1:Domain        % Loop for two domains lying in first coordinate
    % Initialize matrices in x and y axis
    X = zeros(m,n) ;
    Y = zeros(m,n) ;

    for i = 1:m
        Xi = xi(i) ;
        for j = 1:n
            Eta = eta(j) ;
        
            % Transfinite Interpolation 
            XY = (1-Eta)*Xb(Xi,d)+Eta*Xt(Xi,d)+(1-Xi)*Xl(Eta,d)+Xi*Xr(Eta,d)......
                -(Xi*Eta*Xt(1,d)+Xi*(1-Eta)*Xb(1,d)+Eta*(1-Xi)*Xt(0,d)+(1-Xi)*(1-Eta)*Xb(0,d)) ;
    
            X(i,j) = XY(1) ;
            Y(i,j) = XY(2) ;
        
        end
    end
    DX{d} = X ;
    DY{d} = Y ;
end
% Arrange the coordinates for each domain
X1 = DX{1} ; Y1 = DY{1} ;       % Grid for first domain
X2 = DX{2} ; Y2 = DY{2} ;       % Grid for second domain
X = [X1 ;X2(m-1:-1:1,:)] ;      % Merge both the domains
Y = [Y1 ;Y2(m-1:-1:1,:)] ;
% Plot 1/4th of the plate
figure(1)
plotgrid(X,Y) ;
% break
% Plot other domains of plate by imaging coordinates
vec = [1 1 ; -1 1 ; -1 -1 ; 1 -1] ;
figure(2)
for quadrant = 1:4
    plotgrid(vec(quadrant,1)*X,vec(quadrant,2)*Y) ;
    hold on
end
    