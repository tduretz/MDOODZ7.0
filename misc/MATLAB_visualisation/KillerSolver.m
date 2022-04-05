clear, close all, clc
addpath('./_M2Di2_functions/')

path = '/Users/tduretz/REPO_GIT/MDOODZ6.0_PreHolidays/SOURCE/';

SuiteSparse = 1;                                                             % 0=NO SuiteSparse  1=SuiteSparse
noisy       = 3;                                                             % 0=NO Output       1=Std output     2=Very noisy  3=Very noisy + kspgcr output
Newton      = 1; 
cpu = zeros(9,1);
% Non-linear solver settings
nmax_glob   = 200;                                                           % Max. number of non-linear iterations
tol_glob    = 1e-10;                                                          % Global nonlinear tolerance
eps_kspgcr  = 1e-4; 
% Linear solver settings
gamma       = 1/(2.5e4);
nPH         = 30;                                                            % Max. number of linear iterations
tol_linu    = tol_glob/5000;                                                 % Velocity tolerance
tol_linp    = tol_glob/1000;                                                 % Divergence tolerance


filename = [path 'Jacobian_01cpu_step08_iter00.gzip.h5'];
[ J, grad, div, PPI, fu, fp, bu, bp ] = ExtractDecoupledMatrixFromMDoodz( filename );

filename = [path 'Stokes_01cpu_step08_iter00.gzip.h5'];
[ K, ~, ~, ~, ~, ~, ~, ~ ] = ExtractDecoupledMatrixFromMDoodz( filename );


% bu = fu;
% bp = fp;

du = zeros(size(fu));
dp = zeros(size(fp));

fu = -(J*du + grad*dp - bu); 
fp = -(div*du  - bp);

PPI = PPI*(1/gamma);


fprintf('min = %2.7e max= %2.7e\n', min(du(:)), max(du(:)) )
fprintf('min = %2.7e max= %2.7e\n', min(dp(:)), max(dp(:)) )

fprintf('min = %2.7e max= %2.7e\n', min(bu(:)), max(bu(:)) )
fprintf('min = %2.7e max= %2.7e\n', min(bp(:)), max(bp(:)) )

fprintf('min = %2.7e max= %2.7e\n', min(fu(:)), max(fu(:)) )
fprintf('min = %2.7e max= %2.7e\n', min(fp(:)), max(fp(:)) )



if Newton==0, J = K; end
tic
PC  = 1/2*(J'+ J);                                                     % Symmetrisation of Jacobian
% PC  = K;
Jt  =  J - grad*(PPI*div); 
Jts = PC - grad*(PPI*div);                                            % Velocity Schur complement of Jacobian operator
[Jcs,e,s] = chol(Jts,'lower','vector');                               % Choleski factorization of Jts
norm_r = 0; its = 0;
cpu(5)=cpu(5)+toc;
tic
% Powell-Hestenes iterations
fu0 = fu;                                                                    % Save linear norm 0
its_KSP_tot = 0;
for itPH=1:nPH
    fut  = bu - grad*dp - grad*PPI*bp;                                       % Iterative right hand side
    [du,norm_r,its] = M2Di2_kspgcr_m(Jt,fut,du,Jcs,s,eps_kspgcr,noisy,SuiteSparse); % Apply inverse of Schur complement
    its_KSP_tot = its_KSP_tot + its;
    dp   = dp + PPI*(bp -  div*du);                                          % Pressure corrctions
    fprintf('min = %2.7e max= %2.7e\n', min(du(:)), max(du(:)) )
    fprintf('min = %2.7e max= %2.7e\n', min(dp(:)), max(dp(:)) )
    sum(du(:))
    sum(dp(:))

    fu  = bu -   J*du  - grad*dp;                                           % Compute linear velocity residual
    fp  = bp - div*du;                                                      % Compute linear pressure residual
   
    fprintf('min = %2.7e max= %2.7e\n', min(fu(:)), max(fu(:)) )
    fprintf('min = %2.7e max= %2.7e\n', min(fp(:)), max(fp(:)) )

    if noisy>1, fprintf('--- iteration %d --- \n',itPH);
        fprintf('  Res. |u| = %2.2e \n',norm(fu)/length(fu));
        fprintf('  Res. |p| = %2.2e \n',norm(fp)/length(fp));
        fprintf('  KSP GCR: its=%1.4d, Resid=%1.4e \n',its,norm_r); end
    if ((norm(fu)/length(fu)) < tol_linu) && ((norm(fp)/length(fp)) < tol_linp), break; end
    if ((norm(fu)/length(fu)) > (norm(fu0)/length(fu)) && norm(fu)/length(fu) < tol_glob), fprintf(' > Linear residuals do no converge further:\n'); break; end
    fu0 = fu;
end
if noisy>=1, fprintf(' - Linear res. ||res.u||=%2.2e\n', norm(fu)/length(fu) );
    fprintf(' - Linear res. ||res.p||=%2.2e\n', norm(fp)/length(fp) );
    fprintf('  KSP GCR tot: KSP_its_tot=%1.4d \n',its_KSP_tot);
end