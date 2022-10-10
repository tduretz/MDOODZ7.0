%==========================================================================
% Version for Matlab which can use vector multiplication
% Version for C-code which is based on components, because no vector
% multiplication
% Director should not be orthogonal to shear direction
clear all, close all, clc

% Simple shear
dt          = 1e-2;
dudx        = 0;
dudy        = 1;
dvdx        = 0;
dvdy        = 0;

mu_N        = 1;
mu_S        = mu_N/10;
C_ISO       = 2*mu_N*[1 0 0; 0 1 0; 0 0 1];
w12         = 1/2*(dudy-dvdx);
e12         = 1/2*(dudy+dvdx);
L           = [dudx dudy; dvdx dvdy];
W           = [0 w12; -w12 0];
D           = [0 e12; e12 0];
D_vec       = [dudx dvdy e12]';
Tau_vec     = C_ISO*D_vec;
Tauxy       = Tau_vec(3);
Gamma       = 0;
F_inc       = [1 dt*dudy; 0 1];
F           = [1       0; dt*dvdx 1];

Director    = [1; 1];
Norm_dir    = sqrt( Director(1)^2+Director(2)^2 );
Director    = Director ./ Norm_dir;
Director2   = Director;
Norm_dir    = sqrt( Director(1)^2+Director(2)^2 );
Nx          = Director(1);
Ny          = Director(2);
Angle       = atand(Nx/Ny);
Fol_x       = -Ny/Nx;
Fol_y       = 1;

it          = 1;
while Ny(end) > -0.99 %it = 2:500;
    it              = it+1;

    F               = F_inc*F;
    B               = F*F';                 % Left Cauchy-Green tensor
    [BV,BE]         = eigs(B);              % Spectral decomposition
    VE              = diag(sqrt(diag(BE))); % Square roots of eigenvalues
    FS              = BV*VE;                % Scale eigenvectors with eigenvalues
    psa1            = FS(:,2);              % Finite strain axis 1
    psa2            = FS(:,1);              % Finite strain axis 2
    Aratio          = sqrt( psa2(1)^2 + psa2(2)^2 ) / sqrt( psa1(1)^2 + psa1(2)^2 );

    Director_dot    = W*Director - D*Director*(Director'*Director) + (Director'*(D*Director))*Director;
    Director_dot2   = -L'*Director2;

    Director        = Director + Director_dot*dt;
    Director2       = Director2 + Director_dot2*dt;
    Norm_dir(it)    = sqrt(Director(1)^2+Director(2)^2);
    Director        = Director ./ Norm_dir(it);
    a0              = 2*Director(1)^2*Director(2)^2;
    a1              = (Director(1)*Director(2)^3-Director(1)^3*Director(2));
    a1              = Director(1)*Director(2) * (-Director(1)^2 + Director(2)^2);
    C_ANI           = [-a0 a0 2*a1; a0 -a0 -2*a1; a1 -a1 -1+2*a0];

    S_C             = C_ISO + 2*(mu_N-mu_S) * C_ANI;

    Tau_vec         = S_C*D_vec;
    Tauxy(it)       = Tau_vec(3);
    Gamma(it)       = Gamma(it-1) + e12*dt;
    Nx(it)          = Director(1);
    Ny(it)          = Director(2);
    Angle(it)       = atand(Nx/Ny);
    Fol_x(it)       = -Ny(it)/Nx(it);
    Fol_y(it)       = 1;

    if mod(it,5)==0
        subplot(221)
        plot(Gamma,Tauxy,'-ok'), grid on
        xlabel('\gamma')
        ylabel('\tau_{xy}')
        subplot(223)
        plot(Gamma,Nx,'-ok'), hold on,plot(Gamma,Ny,'-sb'),plot(Gamma,Norm_dir,'-+r'), grid on
        legend('N_x','N_y','norm','location','southwest')
        xlabel('\gamma')
        ylabel(' ')

        subplot(222)
        hold off
        plot([0 Nx(it)   ],[0 Ny(it)]   ,'-r','linewidth',2), grid on, hold on
        plot([0 Director2(1)   ]/sqrt(Director2(1)^2+Director2(2)^2),[0 Director2(2)]/sqrt(Director2(1)^2+Director2(2)^2)   ,'-o','linewidth',1)

        plot([Fol_x(it) -Fol_x(it)],[Fol_y(it) -Fol_y(it)],'-b','linewidth',2)
        plot([0 Nx(it)],[0 0],'-r')
        plot([Nx(it) Nx(it)],[0 Ny(it)],'-r')
        quiver(0,0, -psa2(1)*1e2, -psa2(2)*1e2,'color','g','linewidth',1,'AutoScale','off');
        axis equal, axis([-1 1 -1 1])
        xlabel('x')
        ylabel('y')
        subplot(224)
        hold off
        quiver(0,0, psa1(1), psa1(2),'color','b','linewidth',2,'AutoScale','off'); hold on
        quiver(0,0, psa2(1), psa2(2),'color','b','linewidth',2,'AutoScale','off');
        quiver(0,0,-psa1(1),-psa1(2),'color','b','linewidth',2,'AutoScale','off');
        quiver(0,0,-psa2(1),-psa2(2),'color','b','linewidth',2,'AutoScale','off');
        axis equal, axis(2*[-1 1 -1 1])
        title(['Aspect ratio: ',num2str(Aratio,2)])
        set(gcf,'position',[393 97.4000 978 684.4000])
        drawnow
    end
end


% % Expression from Mulhaus
% n1    = Director(1);
% n2    = Director(2);
% ddotx = (-(dudx-dvdy)*n1*n2 - dvdx*n2^2 + dudy*n1^2)*n2;
% ddoty = ( (dudx-dvdy)*n1*n2 + dvdx*n2^2 - dudy*n1^2)*n1;
% ddot  = [ddotx;ddoty]