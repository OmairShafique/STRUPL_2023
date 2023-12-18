function [der,fun] = fmquad(samp,Element_Type,dim, ig,jg,iel)%,lg

%% To be called only for Element_Type == 8 or == 4
%
global ngpb ngps connec geom

% This function
% returns the vector of the shape function and their
% derivatives with respect to xi and eta at the gauss points for
% an 8-nodded quadrilateral
%
% Element_Type=8;
% Element_Type=3;
% Element_Type=4;
% dim=2;
% dim=3;
xi=samp(ig,1);
eta=samp(jg,1);
% zeta=samp(lg,1);
etam=(1.-eta);
etap=(1.+eta);
xim=(1.-xi);
xip=(1.+xi);
% zetam=(1.-zeta);
% zetap=(1.+zeta);

%

switch (Element_Type)

    case (3) % for triangular element
        if Element_Type==3 && dim==2 && ngpb ~=0 && ngps ~=0
            % shape functions
            fun(1) = 1. - xi - eta;
            fun(2) =  xi;
            fun(3) =  eta;
            % derivatives


            der(1,1)= -1.; der(1,2)=-1;
            der(2,1)= 1.; der(2,2)= 0;
            der(3,1)=0; der(3,2)=1;
        else
            x1 = geom(connec(iel,1),1); y1 = geom(connec(iel,1),2);
            x2 = geom(connec(iel,2),1); y2 = geom(connec(iel,2),2);
            x3 = geom(connec(iel,3),1); y3 = geom(connec(iel,3),2);

            A = (0.5)*det([1 x1 y1; ...
                1 x2 y2; ...
                1 x3 y3]);

            m11 = (x2*y3 - x3*y2)/(2*A); % reanme to ders
            m21 = (x3*y1 - x1*y3)/(2*A);
            m31 = (x1*y2 - y1*x2)/(2*A);
            m12 = (y2 - y3)/(2*A);
            m22 = (y3 - y1)/(2*A);
            m32 = (y1 - y2)/(2*A);
            m13 = (x3 - x2)/(2*A);
            m23 = (x1 - x3)/(2*A);
            m33 = (x2 -x1)/(2*A);

            % der(1,1)  = m11;
            % der(2,1)  = m21;
            % der(3,1)  = m31;
            der(1,1)  = m12;
            der(1,2)  = m22;
            der(2,1)  = m32;
            der(2,2)  = m13;
            der(3,1)  = m23;
            der(3,2)  = m33;

            fun = [(m11+m12+m13)  0;...
                0  (m11+m12+m13);...
                (m21+m22+m23)  0;...
                0  (m21+m22+m23);...
                (m31+m32+m33)  0;...
                0  (m31+m32+m33)];

        end
    case (4) % for rectangular element

        % shape functions
        fun(1)=0.25*xim*etam;
        fun(2)=0.25*xip*etam;
        fun(3)=0.25*xip*etap;
        fun(4)=0.25*xim*etap;
        % derivatives
        der(1,1)=-0.25*etam; der(1,2)=-0.25*xim;
        der(2,1)=0.25*etam;  der(2,2)=-0.25*xip;
        der(3,1)=0.25*etap;  der(3,2)=0.25*xip;
        der(4,1)=-0.25*etap; der(4,2)=0.25*xim;

    case (8) % for 8-noded Element

        if Element_Type==8 && dim==2
            % shape functions
            fun(1) = -0.25*xim*etam*(1.+ xi + eta);
            fun(2) = 0.5*(1.- xi^2)*etam;
            fun(3) = -0.25*xip*etam*(1. - xi + eta);
            fun(4) = 0.5*xip*(1. - eta^2);
            fun(5) = -0.25*xip*etap*(1. - xi - eta);
            fun(6) = 0.5*(1. - xi^2)*etap;
            fun(7) = -0.25*xim*etap*(1. + xi - eta);
            fun(8) = 0.5*xim*(1. - eta^2);
            % derivatives
            der(1,1)=0.25*etam*(2.*xi + eta); der(1,2)=-1.*etam*xi;
            der(1,3)=0.25*etam*(2.*xi-eta); der(1,4)=0.5*(1-eta^2);
            der(1,5)=0.25*etap*(2.*xi+eta); der(1,6)=-1.*etap*xi;
            der(1,7)=0.25*etap*(2.*xi-eta); der(1,8)=-0.5*(1.-eta^2);
            %
            der(2,1)=0.25*xim*(2.*eta+xi); der(2,2)=-0.5*(1. - xi^2);
            der(2,3)=-0.25*xip*(xi-2.*eta); der(2,4)=-1.*xip*eta;
            der(2,5)=0.25*xip*(xi+2.*eta); der(2,6)=0.5*(1.-xi^2);
            der(2,7)=-0.25*xim*(xi-2.*eta); der(2,8)=-1.*xim*eta;
        else
            if Element_Type==8 && dim==3  %brick element
                % shape functions
                fun(1) = 0.125*xim*etam*zetam;
                fun(2) = 0.125*xip*etam*zetam;
                fun(3) = 0.125*xip*etap*zetam;
                fun(4) = 0.125*xim*etap*zetam;
                fun(5) = 0.125*xim*etam*zetap;
                fun(6) = 0.125*xip*etam*zetap;
                fun(7) = 0.125*xip*etap*zetap;
                fun(8) = 0.125*xim*etap*zetap;
                % derivatives
                der(1,1)=0.125*etam*zetam; der(1,2)=0.125*etam*zetam;
                der(1,3)=0.125*etap*zetam; der(1,4)=0.125*etap*zetam;
                der(1,5)=0.125*etam*zetap; der(1,6)=0.125*etam*zetap;
                der(1,7)=0.125*etap*zetap; der(1,8)=0.125*etap*zetap;
                %
                der(2,1)=0.125*xim*zetam; der(2,2)=0.125*xip*zetam;
                der(2,3)=0.125*xip*zetam; der(2,4)=0.125*xim*zetam;
                der(2,5)=0.125*xim*zetap; der(2,6)=0.125*xip*zetap;
                der(2,7)=0.125*xip*zetap; der(2,8)=0.125*xim*zetap;
                %
                der(3,1)=0.125*xim; der(3,2)=0.125*xip;
                der(3,3)=0.125*xip; der(3,4)=0.125*xim;
                der(3,5)=0.125*xim; der(3,6)=0.125*xip;
                der(3,7)=0.125*xip; der(3,8)=0.125*xim;
            end
        end
end
end