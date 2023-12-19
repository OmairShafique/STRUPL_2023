% THIS PROGRAM USES AN 3-NODE LINEAR TRIANGULAR ELEMENT FOR THE
% LINEAR ELASTIC STATIC ANALYSIS OF A TWO DIMENSIONAL PROBLEM
%
clear all
clc
tic
%
% Make these variables global so they can be shared by other functions
%
global nnd nel nne nodof eldof n
global geom connec dee nf Nodal_loads
%
format long g
%
% ALTER NEXT LINES TO CHOOSE THE NAME OF THE OUTPUT FILE
%
fid =fopen('CST_COARSE_MESH_RESULTS.txt','w');
%
% To change the size of the problem or change elastic properties
% supply another input file
%
CST_COARSE_MESH_DATA;
PlotMesh(geom,connec)
%
%%%%%%%%%%%%%%%%%%%%%%%%%% End of input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Assemble the global force vector
%
fg=zeros(n,1);
for i=1: nnd
if nf(i,1) ~= 0
fg(nf(i,1))= Nodal_loads(i,1);
end
if nf(i,2) ~= 0
fg(nf(i,2))= Nodal_loads(i,2);
end
end
%
% Assembly of the global stiffness matrix
%
% initialize the global stiffness matrix to zero
%
KK = zeros(n, n);
%
for i=1:nel
[bee,g,A] = elem_T3(i); % Form strain matrix, and steering vector
ke=thick*A*bee'*dee*bee; % Compute stiffness matrix
KK=form_KK(KK,ke, g); % assemble global stiffness matrix
end
%
%
%%%%%%%%%%%%%%%%%%%%%%% End of assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Compute global displacement vector and reaction%%
[delta,Reaction,node_displacement,SIGMA,STRAIN]=Elastic_solve(KK,fg,nf,dee,bee,thick);
%
type='average';
type='sum';
stressavg = zeros(nnd,3);
for ind = 1:nnd
    [r c] = find(connec==ind) ;
    switch type
        case 'average'
        share = length(r) ;
        case 'sum'
        share = 1 ;
    end
    ele = 4*(r-1)+c ;
    stress=[sum(stress_N(ele,1))/share.....
            sum(stress_N(ele,2))/share.....
            sum(stress_N(ele,3))/share];                        
    stressavg(ind,:)=stress;
end

% %
%% Compute for each element the plastic multiplier that activates the plastic mode This part is a function related to yield function drucker pragher [alpha_b,Normal_stress]
%Macro-cracks parameters
%
sigma_t = 0.1198; % Limit tension [MPa]
tau = 0.1198; % Limit shear [MPa]
tau_m = tau*3 ; %Valore di tau per macrofess
Gf_t = 0.035; % Fracture Energy traction [N/mm]
Gf_s = Gf_t ; % Fracture Energy shear[N/mm]
w_t = 2 * Gf_t/sigma_t ;
w_s = 2 * Gf_s/(tau) ;
h_s_t = -sigma_t /w_t  ; %softening parameter
h_s_s = -tau /w_s  ; %softening parameter
%Micro-cracks (Drucker-Prager parameters)
%
Plastic_mode = [1:1:nel]'; %Plastic mode
d_u = 0.019 ; %Ultimate deformation
coeff_h = 0.01 ; %percentuale di resistanza in corrispondenza del quale si attiva l'hardning
% h_t = ((1-coeff_h)*sigma_t /d_u) ; %hardening parameter
% h_s = ((1-coeff_h)*tau/d_u) ; %hardening parameter
h_t = (0.9*sigma_t /d_u) ; %hardening parameter
sigma_t_r = coeff_h * sigma_t ; % Limit traction
sigma_c = sigma_t_r  ; % Limit compression
%-------------------------
% sigma_t_tol = sigma_t_r + eps*sigma_t_r ;
% sigma_c_tol = sigma_c + eps*sigma_c ;
Aa = (2/sqrt(3))*((sigma_t_r*sigma_c)/(sigma_t_r+sigma_c)) ;
Bb = (1/sqrt(3))*((sigma_t_r-sigma_c)/(sigma_t_r+sigma_c)) ;
% Aa_tol = (2/sqrt(3))*((sigma_t_tol*sigma_c_tol)/(sigma_t_tol+sigma_c_tol)) ;
% Bb_tol = (1/sqrt(3))*((sigma_t_tol-sigma_c_tol)/(sigma_t_tol+sigma_c_tol)) ;

for j = 1 : size(Plastic_mode,1)
    i = Plastic_mode(j) ;
    %Compute the principal stress 
    SIGMA_principal(i,1) = ((SIGMA(i,1)+SIGMA(i,2))/2) + (((SIGMA(i,1)-SIGMA(i,2))/2)^2 + (SIGMA(i,3)^2))^(1/2) ; %sigmaI
    SIGMA_principal(i,2) = ((SIGMA(i,1)+SIGMA(i,2))/2) - (((SIGMA(i,1)-SIGMA(i,2))/2)^2 + (SIGMA(i,3)^2))^(1/2) ; %sigmaII
    SIGMA_principal(i,3) = (((SIGMA(i,1)-SIGMA(i,2))/2)^2 + (SIGMA(i,3)^2))^(1/2) ; %taumax
    
    %Compute the first stress invariant
    I1(j,1) = SIGMA(i,1) + SIGMA(i,2) ;
    
    %Compute the  second invariant of the deviatoric stress it describes a state of pure shear.
    J2(j,1) = 1/6*(SIGMA(i,1)-SIGMA(i,2))^2 + 1/6*SIGMA(i,1)^2 + 1/6*SIGMA(i,2)^2 + SIGMA(i,3)^2 ;
    alpha_b(j,1) = Aa/(sqrt(J2(j))-Bb*I1(j)) ;
end
clear i j 

    % Solve the elastic problem amplifying the load  by multiplication of
    % fg *alpha to computes the normal stress
    %     [u_alpha]=solve(Kk,alpha_b(j)*F,bc) ;
   

eld_alpha = []; %element stress matrix
Normal_stress = [] ;
SIGMA_alpha = [] ;
for j = 1 : size(Plastic_mode,1)
    delta_alpha = KK\fg*alpha_b(j) ; % solve for unknown displacements
    Reaction_alpha=KK*delta_alpha-fg*alpha_b(j);
    %
    node_disp_alpha=zeros(nnd,2);
    %
    for i=1: nnd %
        if nf(i,1) == 0 %
            x_disp =0.; %
        else
            x_disp = delta_alpha(nf(i,1)); %
        end
        %
        if nf(i,2) == 0 %
            y_disp = 0.; %
        else
            y_disp = delta_alpha(nf(i,2)); %
        end
        node_disp(i,:) =[x_disp y_disp];
    end
    %
    % Retrieve the x_coord and y_disp of the nodes located on the neutral axis
    %
    k = 0;
    for i=1:nnd
        if geom(i,2)== 0.
            k=k+1;
            x_coord(k) = geom(i,1);
            vertical_disp(k)=node_disp(i,2);
        end
    end
    %% Compute stress and strain%% 
    %
    for i=1:nel
        [bee,g,A] = elem_T3(i); % Form strain matrix, and steering vector
        eld_alpha=zeros(eldof,1); % Initialize element displacement to zero
        for m=1:eldof
            if g(m)==0
                eld_alpha(m)=0.;
            else %
                eld_alpha(m)=delta_alpha(g(m)); % Retrieve element displacement
            end
        end
        %
        eps_alpha=bee*eld_alpha; % Compute strains
        EPS_alpha(i,:)=eps_alpha ; % Store strains for all elements
        sigma_alpha=dee*eps_alpha; % Compute stresses
        SIGMA_alpha(i,:)=sigma_alpha ; % Store stress for all elements
    end
    clear i

    i = Plastic_mode(j) ;

    %Compute the  second invariant of the deviatoric stress it describes a state of pure shear.
    J2_alpha = 1/6*(SIGMA_alpha(i,1)-SIGMA_alpha(i,2))^2 + 1/6*SIGMA_alpha(i,1)^2 + 1/6*SIGMA_alpha(i,2)^2 + SIGMA_alpha(i,3)^2 ;

    %Compute normal at a surface point 
    dphi_dsigmax = -Bb + 0.5*J2_alpha^(-0.5) * (1/3*(SIGMA_alpha(i,1)-SIGMA_alpha(i,2))+1/3*SIGMA_alpha(i,1)) ;
    dphi_dsigmay = -Bb + 0.5*J2_alpha^(-0.5) * (-1/3*(SIGMA_alpha(i,1)-SIGMA_alpha(i,2))+1/3*SIGMA_alpha(i,2)) ;
    dphi_dsigmaxy =  0.5*J2_alpha^(-0.5) * 2*SIGMA_alpha(i,3) ;
    %
    nnn = [dphi_dsigmax dphi_dsigmay dphi_dsigmaxy]  ;
    nnn_n = nnn/norm(nnn) ;
    Normal_stress = [Normal_stress; nnn_n] ;
%Implement the recovery  stress to find the Normal_nodal_stress
% Avverage for nodes
% Normal_stressR=zeros(nnd,6); % Initial global stress matrix by nodes
% for i=1:nnd;
%  Normal_stressR(i,1:3)=zeros(3,1);
%  count=0;
%  for j=1:nel;
%  if ismember(i,ET3(j,:));
%  count=count+1;
%  Normal_stressR(i,1:3)=Normal_stressR(i,1:3)+ Normal_stress; % get sum
%  end;
%  end;
%  Normal_stressR(i,1:3)=Normal_stressR(i,1:3)/count; % average
%  Normal_stressR(i,4)=0.5*(Normal_stressR(i,1)+Normal_stressR(i,2))+0.5*sqrt((Normal_stressR(i,1)-Normal_stressR(i,2))^2+4*Normal_stressR(i,3)^2); % max main normal stress
%  disp(StressR(i,4));
%  Normal_stressR(i,5)=0.5*(Normal_stressR(i,1)+StressR(i,2))-0.5*sqrt((Normal_stressR(i,1)-Normal_stressR(i,2))^2+4*Normal_stressR(i,3)^2); % min main normal stress
%  Normal_stressR(i,6)=0.5*(Normal_stressR(i,4)-Normal_stressR(i,5)); % max tangential stress
% end;
% end

end
clear j
%%  Computes Resistance vector Function%% 
crack_path_2=[];
% this vector is the former IL vector from manuela and franchi codes is formed by 1st column prone crack nodes and 2dn column is the length 
% a pre-check loop is replaced here where it  checks the limit tension if
% exceeds in the SIGMA_alpha array.
% crack_path=[1 5;1 5; 2 10; 2 10];
% working on this part that is crucial 
crack_path_2 = Gcrack(Normal_stress, nnd, connec, sigma_t, geom,'Crack_path');
idx=crack_path_2(:,3)~=0; 
Number_of_crack_nodes=sum(idx(:));
crack_direction=repelem(crack_path_2(:,3),2,1);
crack_nodes_lines=[ 2 5 7 0 0 0 0;2 5  8 10 0 0 0;2 5 8 11 13 0 0;2 5 8 11 14 16 0;2 5 8 11 14 17 0 ;2 5 8 11 14 17 19;3 5 7 0 0 0 0;3 5 8 10 0 0 0;3 5 8 11 13 0 0;3 5 8 11 14 16 0;3 5 8 11 14 17 0;3 5 8 11 14 17 19];
Number_of_crack_lines= length(crack_nodes_lines);
%
% The volume of elements is important in the integration
Volume = [];
for i = 1 : Plastic_mode
    Volume = [Volume; polyarea(ex(Plastic_mode(i),:),ey(Plastic_mode(i),:))*thick] ;
end
RR = zeros(size(crack_direction));
RR(1:2:end) = crack_direction(1:2:end)*thick*sigma_t;% find the resistance of nodes that failure in tension
RR(2:2:end) = crack_direction(2:2:end)*thick*tau_m;%nodes that failure in shear

% for crack nodes the resistance is the multplication of stress by normal
% resistance relative ai modi di fessurazione localizzata (prodotto della
% % tensione  limite per l'area di influenza di ciascun nodo, calcolata a

Resistance_cracking = zeros(2*2*size(crack_path_2(:,1),1),1) ;
for rr = 1 : size(RR,1)
    Resistance_cracking(2*rr-1) = RR(rr) ;
end
clear rr

Resistance_plastic = [] ;
for j = 1 : size(Plastic_mode,1)
    i = Plastic_mode(j) ;
    RR_p = SIGMA_alpha(j,:)* Normal_stress(j,:)' ;
    Resistance_plastic = [Resistance_plastic; RR_p] ;
end
clear j
%Moltiplico per il volume di ciscun elemento (intragrazione nel volume)

Resistance_plastic = Resistance_plastic.*Volume ;

Resistance_vector = [Resistance_cracking; Resistance_plastic];
%% Assemble vector b %%%%%%%%
%Assemblage of nodal rotation matrix
% beta
% for cc = 1 : size(Cp,1)
%     if cc == 1
%         [NN] = Rotation_matrix(Beta(cc,2),Beta(cc,3)) ;
%     else
%         [NN] = [NN; Rotation_matrix(Beta(cc,2),Beta(cc,3))] ;
%     end
% end
% Computes internal force
Vector_b_plastic= [];
Vector_b_cracking= [];
for j = 1 : size(Plastic_mode,1)
     SIGMA_principal = SIGMA(Plastic_mode(j),:) * Normal_stress(Plastic_mode(j),:)' ;
    Vector_b_plastic = [Vector_b_plastic; SIGMA_principal] ;
end
clear j
%Moltiplico per il volume di ciscun elemento (intragrazione nel volume)
Vector_b_plastic = Vector_b_plastic.*Volume; 

% La componente associata alla fessurazione localizzata bb_f si calcola con
% apposita funzione e viene poi modificata per tenere conto del softening

% [Vector_B] = vector_b(et,ex,ey,N,C,F,b,bc,D,Cp,Cpp,Be,Bee,NN) ;
% 
% bb_f = zeros(2*size(bbf,1),1) ;
% for xx = 1 : size(bbf,1)
%     bb_f(2*xx-1) = bbf(xx);
%     bb_f(2*xx) = bbf(xx);
% end
% clear xx
% 
% Vector_B = [bb_f; Vector_b_plastic] ;
%

%% Assemble Matrix A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% [A_Matrix]=form_Matrix_A(Normality_Matrix,dee,bee,KK);
% Print results to file
%
print_CST_results;
%
% Plot the stresses in the x_direction
%
x_stress = SIGMA(:,1);
cmin = min(x_stress);
cmax = max(x_stress);
clim([cmin cmax])
patch('Faces', connec, 'Vertices', geom, 'FaceVertexCData',x_stress, ...
'Facecolor','flat','Marker','.')
colorbar;
% % Deformed --------------------------------------------------------------
% figure('Name','Deformed shape')
% k=2;
% x=ex' ;
% y=ey';
% xc=[x; x(1,:)];
% yc=[y; y(1,:)];
% x_d=(ex+k*ed(:,[1 3 5]))';   
% y_d=(ey+k*ed(:,[2 4 6]))';
% xc_d=[x_d; x_d(1,:)];
% yc_d=[y_d; y_d(1,:)];
% hold on
% axis equal
% plot(xc,yc,'Color',[0.7 0.7 0.7]);
% plot(xc_d,yc_d,'k') ;
% set(0,'defaultlinelinewidth',1.5);
% ax.GridLineStyle = '-';
% % grid on
% box on
% % xlim([-200 1400]) ;
% % ylim([-200 1400]) ;
% saveas(gcf,'figure18.png')
toc