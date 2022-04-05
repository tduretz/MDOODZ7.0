% function Main

clear 
close all
clc

% addpath '/Users/tduretz/SuiteSparse4.4.3/MATLAB_Tools/'
% addpath '/Users/tduretz/SuiteSparse4.4.3/CHOLMOD/MATLAB/'
% addpath '/Users/tduretz/SuiteSparse4.4.3/CXSparse/MATLAB/CSparse/'
% addpath '/Users/tduretz/SuiteSparse4.4.3/UMFPACK/MATLAB/'
% addpath '/Users/tduretz/SuiteSparse4.4.3/AMD/MATLAB/'
% addpath '/Users/tduretz/SuiteSparse4.4.3/CAMD/MATLAB/'
% addpath '/Users/tduretz/SuiteSparse4.4.3/CCOLAMD/MATLAB/'

% here you change this parameter to 16, 32 or 64 (resolution)
res = 21;

solve = 0;

istep = 1;

tags = 0;
visc = 0;
topo = 0;

jac = 0;

path = '/Users/tduretz/REPO_GIT/MDOODZ6.0/SOURCE/';


filename = [path 'MatrixCoupled.gzip_1cpu.h5'];
% filename = [path 'Matrix.gzip_1cpu.h5'];

Ic  = hdf5read(filename,'/matrix/I');
J   = hdf5read(filename,'/matrix/J');
Ic  = Ic+1;
J   = J+1;
V1  = hdf5read(filename,'/matrix/V');
L   = hdf5read(filename,'/matrix/rhs');

if jac == 1
    filename_jac = [path 'Jacobian.gzip_4cpu.h5'];
    IcJ  = hdf5read(filename_jac ,'/matrix/I');
    JJ   = hdf5read(filename_jac ,'/matrix/J');
    IcJ  = IcJ+1;
    JJ   = JJ+1;
    V1J  = hdf5read(filename_jac ,'/matrix/V');
    LJ   = hdf5read(filename_jac ,'/matrix/rhs');
end

equ = hdf5read(filename,'/numbering/eqn_u');
eqv = hdf5read(filename,'/numbering/eqn_v');
eqp = hdf5read(filename,'/numbering/eqn_p');

tag_s = hdf5read(filename,'/fields/tag_s');
tag_n = hdf5read(filename,'/fields/tag_n');
tag_u = hdf5read(filename,'/fields/tag_u');
tag_v = hdf5read(filename,'/fields/tag_v');

eta_s = hdf5read(filename,'/fields/eta_s');
eta_n = hdf5read(filename,'/fields/eta_n');

params = hdf5read(filename,'/model/params');

nx = params(1);
nz = params(2);
dx = params(3);
dz = params(4);

% Read mesh
filename0 = [path,'Output',num2str(0,'%05d'),'','.gzip.h5'];
% filename0 = [path,'DEBUG',num2str(istep,'%05d'),'','.gzip.h5'];
params     = hdf5read(filename0,'/Model/Params');
xg_coord   = hdf5read(filename0,'/Model/xg_coord');  xg_coord  = cast(xg_coord, 'double');
zg_coord   = hdf5read(filename0,'/Model/zg_coord');  zg_coord  = cast(zg_coord, 'double');
xc_coord   = hdf5read(filename0,'/Model/xc_coord');  xc_coord  = cast(xc_coord, 'double');
zc_coord   = hdf5read(filename0,'/Model/zc_coord');  zc_coord  = cast(zc_coord, 'double');
xvz_coord  = hdf5read(filename0,'/Model/xvz_coord'); xvz_coord = cast(xvz_coord, 'double');
zvx_coord  = hdf5read(filename0,'/Model/zvx_coord'); zvx_coord = cast(zvx_coord, 'double');

% Read chain
filename_p = [path,'Particles',num2str(istep,'%05d'),'','.gzip.h5'];

% filename_p = [path,'DEBUGParticles',num2str(istep,'%05d'),'','.gzip.h5'];

if topo
    xtopo = hdf5read(filename_p,'/Topo/x');
    ztopo = hdf5read(filename_p,'/Topo/z');
    vxtopo = hdf5read(filename_p,'/Topo/vx');
    vztopo = hdf5read(filename_p,'/Topo/vz');
    xtopo = cast(xtopo, 'double');
    ztopo = cast(ztopo, 'double');
    height = hdf5read(filename_p,'/Topo/height');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompress index I
I = zeros(size(J));

m = 0;

I(1) = 1;
p    = 1;

for k=2:length(Ic)
    
    off = Ic(k) - Ic(k-1);
    
    for l=1:off
        I(m+l) = p;
    end
    m = m + off;
    p = p + 1;
 
end

rhs = L;
I1  = cast(I, 'double');
J1  = cast(J, 'double');
M    = sparse(I1, J1, V1);

if jac == 1
    
    % Decompress index I
    IJ = zeros(size(JJ));
    
    m = 0;
    
    IJ(1) = 1;
    p    = 1;
    
    for k=2:length(IcJ)
        
        off = IcJ(k) - IcJ(k-1);
        
        for l=1:off
            IJ(m+l) = p;
        end
        m = m + off;
        p = p + 1;
        
    end
    
    rhsJ = LJ;
    I1J  = cast(IJ, 'double');
    J1J  = cast(JJ, 'double');
    MJ    = sparse(I1J, J1J, V1J);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Read Jacobian of Matlab code
% Jmat = load('/Users/tduretz/Desktop/NEWTON_STOKES/MatlabJacobian');

figure(10)
% subplot(131)
% spy(M2'-M2);
% title('M2^t-M2')
% subplot(132)
subplot(121)
spy(M);
title('M')
subplot(122)
% spy(M-M');
title('M-M'' ')
% subplot(122)
% spy(MJ);
% title('M2')
% subplot(133)
% spy(Jmat.J);
% title('M2_{matlab}')

equ = equ(:);
eqv = eqv(:);
eqp = eqp(:);
ivxmin = min(equ(equ~=-1))+1;
iprmin = min(eqp(eqp~=-1))+1;
iprmax = max(eqp(eqp~=-1))+1;

save( 'CoupledMatrix' , 'M', 'L')


if solve == 1
    
    disp('********** Solving Stokes system **********')
    
    disp('1. backslash')
    
    tic 
    x=mldivide(M,L);
    toc    
       min(x)
    max(x)

    tic
    [x, Info] = umfpack2 (M, '\', L) ;
    toc
    Info
    
    min(x)
    max(x)
    
%     spparms('spumoni',2);
    tic
    x=M\L;
    toc
       min(x)
    max(x)
%     
%     % 2. symmetrised matrix backslash
%     %     tic
%     %     x=(M'*M)\(M'*L);
%     %     toc
%     %     min(x)
%     %     max(x)
%      
      
    disp('3. Powell-Hestenes')
    
    Msym=M;
    Lsym=L;
    rSt  = Lsym(1:iprmin-1);
    rPt  = Lsym(iprmin:end);
    a=Msym(ivxmin:iprmin-1,ivxmin:iprmin-1);
    b=Msym(ivxmin:iprmin-1,iprmin:iprmax);
    e=Msym(iprmin:iprmax,ivxmin:iprmin-1);
    f=Msym(iprmin:iprmax,iprmin:iprmax);
    
    tic
    nDirP = 4;
    nit = 3;
    Nx = nx;
    Ny = nz;
    
    tic
    % Assemble weakly compressible preconditionner
    gamma  = 1e7;
    f1inv  = -gamma*speye((Ny-1)*(Nx-1)-nDirP);
    
    % Powell-Hestenes iterations
    nit = 3;
    p   = zeros((Ny-1)*(Nx-1)-nDirP, 1);
    u   = zeros((Ny-1)*(Nx-2) + (Ny-2)*(Nx-1) , 1);
    a1  = a - b*(f1inv*e);
    
    % SuiteSparse   
    tic
    perm = amd2(a1);
    a1   = cs_transpose(a1);
    a1   = cs_symperm(a1,perm);
    a1   = cs_transpose(a1);
    L1   = lchol(a1);
    
    for it=1:nit
        rSt1 = rSt - b*p - b*f1inv*rPt;
        u(perm) = cs_ltsolve( L1, cs_lsolve(L1,rSt1(perm) ) );
        dp   = f1inv*(rPt-e*u);
        p    = p + dp;
        meanDiv = mean(e*u-rPt);
        disp(['it. ' num2str(it) ' Divergence = ', num2str(norm((e*u-rPt)-meanDiv)/length(p) )]);
    end
    
    time2 = toc;
    fu = a*u + b*p - rSt;
    fp = e*u + f*p - rPt;
    disp(['Res. |fu| = ', num2str(norm(fu)/length(fu))]);
    disp(['Res. |fp| = ', num2str(norm(fp)/length(fp))]);
    disp(['Time elapsed: ', num2str(time2), ' s'])
    x = [u; p];
    
%     min(x)
%     max(x)
%     
%     
%     figure(10)
%     subplot(221)
%     spy(a)
%     subplot(222)
%     spy(b)
%     subplot(223)
%     spy(e)
%     subplot(224)
%     spy(f)
    
    if jac == 1
        
        
        disp('********** Solving Jacobian system **********')
        %%
        disp('1. backslash')
        tic
%          spparms('spumoni',2);
        x=MJ\LJ;
        time2 = toc;
        disp(['Time elapsed: ', num2str(time2), ' s'])
        min(x)
        max(x)
        
%         %% 2. symmetrised matrix backslash
% %         tic
% %         x=(MJ'*MJ)\(MJ'*LJ);
% %         time2 = toc;
% %         disp(['Time elapsed: ', num2str(time2), ' s'])
% %         min(x)
% %         max(x)
%           
%         %%
%         disp('3. Powell-Hestenes')
        MsymJ=MJ;
        LsymJ=LJ;
        rStJ  = LsymJ(1:iprmin-1);
        rPtJ  = LsymJ(iprmin:end);
        aJ=MsymJ(ivxmin:iprmin-1,ivxmin:iprmin-1);
        bJ=MsymJ(ivxmin:iprmin-1,iprmin:iprmax);
        eJ=MsymJ(iprmin:iprmax,ivxmin:iprmin-1);
        fJ=MsymJ(iprmin:iprmax,iprmin:iprmax);
%         
%         nDirP = 4;
%         Nx    = nx;
%         Ny    = nz;
%        
%         % Assemble weakly compressible preconditionner
%         gamma  = 1e7;
%         f1inv  = -gamma*speye((Ny-1)*(Nx-1)-nDirP);
% 
%         % Powell-Hestenes iterations
%         nit = 3;
%         p   = zeros((Ny-1)*(Nx-1)-nDirP, 1);
%         u   = zeros((Ny-1)*(Nx-2) + (Ny-2)*(Nx-1) , 1);
%         tic
%         a1J1  = aJ - bJ*(f1inv*eJ);
% 
%         for it=1:nit
%             rSt1J = rStJ - bJ*p - bJ*f1inv*rPtJ;
%             u    = a1J1\rSt1J; 
%             dp   = f1inv*(rPtJ-eJ*u);
%             p    = p + dp;
%             meanDiv = mean(eJ*u);
%             disp(['it. ' num2str(it) ' Divergence = ', num2str(norm((eJ*u)-meanDiv)/length(p) )]);
%         end
%         
%         time2 = toc;
%         fu = aJ*u + bJ*p - rStJ;
%         fp = eJ*u + fJ*p - rPtJ;
%         disp(['Res. |fu| = ', num2str(norm(fu)/length(fu))]);
%         disp(['Res. |fp| = ', num2str(norm(fp)/length(fp))]);
%         disp(['Time elapsed: ', num2str(time2), ' s'])
%         x = [u; p];
% %         min(x)
% %         max(x)
%         
%         %%
        disp('4. Richardson iterations')
        
        tic
        du  = zeros((Ny-1)*(Nx-2) + (Ny-2)*(Nx-1) , 1);
        u   = zeros((Ny-1)*(Nx-2) + (Ny-2)*(Nx-1) , 1);
        p   = zeros((Ny-1)*(Nx-1)-nDirP, 1);
        a1J2 = (aJ) - bJ*(f1inv*eJ);
       
        x      = [u;p];        
        b_zero =  0*speye( double(iprmin-1), (Ny-1)*(Nx-1)-nDirP);
        e_zero =  0*speye((Ny-1)*(Nx-1)-nDirP, double(iprmin-1));
        MsymJ2 = [aJ bJ; eJ, fJ ]; % full Jacobian
%         MsymJ2 = [aJsym, b; e, f1inv ];
                
        F = LsymJ - MsymJ*x;

        for it=1:1000
            dx = M\F;
            x = x + 1.0*dx;
            F = LsymJ - MsymJ*x;
            disp(['Res. |F| = ', num2str(norm(F)/length(F))]);
        end
%         min(x)
%         max(x)
             
        time2 = toc;
        disp(['Time elapsed: ', num2str(time2), ' s'])
        
        %%
        disp('5. Powell-Hestenes')
        nDirP = 4;
        Nx    = nx;
        Ny    = nz;
        
        MsymJ = MJ;
        LsymJ = LJ;
        rStJ  = LsymJ(1:iprmin-1);
        rPtJ  = LsymJ(iprmin:end);
        
        aJ=MsymJ(ivxmin:iprmin-1,ivxmin:iprmin-1);
        bJ=MsymJ(ivxmin:iprmin-1,iprmin:iprmax);
        eJ=MsymJ(iprmin:iprmax,ivxmin:iprmin-1);
        fJ=MsymJ(iprmin:iprmax,iprmin:iprmax);
        
        onesmat = a>0;
        aJsym = onesmat.*aJ;
       
        % Assemble weakly compressible preconditionner
        gamma  = 1e7;
        f1inv  = -gamma*speye((Ny-1)*(Nx-1)-nDirP);

        % Powell-Hestenes iterations
        nit  = 3;
        p    = zeros((Ny-1)*(Nx-1)-nDirP, 1);
        u    = zeros((Ny-1)*(Nx-2) + (Ny-2)*(Nx-1) , 1);
        
        tic
        % Schur complement velocity
%         a1J1 = aJ - bJ*(f1inv*eJ);
          a1J1 = a - bJ*(f1inv*eJ);

        
        % SuiteSparse
        u   = zeros((Ny-1)*(Nx-2) + (Ny-2)*(Nx-1) , 1);
        p   = zeros((Ny-1)*(Nx-1)-nDirP, 1);
        
        tic
%         a1J1 = a1J1'*a1J1;
%         perm = amd2(a1J1);
%         a1J1 = cs_transpose(a1J1);
%         a1J1 = cs_permute(a1J1,perm,perm);
%         a1J1 = cs_transpose(a1J1);
        
        perm = 1:size(rStJ);
        fprintf( 'Reordering took: %2.4f s \n', toc)
%         
%         tic
%         L1 = lchol (a1J1);
%         fprintf( 'UMFPACK took: %2.4f s \n', toc)
% 
%         tic
%         for it=1:nit
%             rSt1J = rStJ - bJ*p - bJ*(f1inv*rPtJ);
%             
%               u(perm) = cs_ltsolve( L1, cs_lsolve(L1,a1J1'*rSt1(perm) ) );
% 
% %             u   = cs_usolve( U1, cs_lsolve(L1, P1*rSt1J(perm))  );
% %             u(perm)    = Q1*u;
%             
%             dp   = f1inv*(rPtJ-eJ*u);
%             p    = p + dp;
%             meanDiv = mean(eJ*u);
%             disp(['it. ' num2str(it) ' Divergence = ', num2str(norm((eJ*u)-meanDiv)/length(p) )]);
%         end
%         
%         perm = 1:size(rStJ);
%         fprintf( 'Reordering took: %2.4f s \n', toc)
%         
        tic
        [L1, U1, P1, Q1] = umfpack2 (a1J1);
        fprintf( 'UMFPACK took: %2.4f s \n', toc)

        tic
        for it=1:nit
            rSt1J = rStJ - bJ*p - bJ*(f1inv*rPtJ);
            
            u   = cs_usolve( U1, cs_lsolve(L1, P1*rSt1J(perm))  );
            u(perm)    = Q1*u;
            
            dp   = f1inv*(rPtJ-eJ*u);
            p    = p + dp;
            meanDiv = mean(eJ*u);
            disp(['it. ' num2str(it) ' Divergence = ', num2str(norm((eJ*u)-meanDiv)/length(p) )]);
        end
        
        fu = aJ*u + bJ*p - rStJ;
        fp = eJ*u + fJ*p - rPtJ;
        disp(['Res. |fu| = ', num2str(norm(fu)/length(fu))]);
        disp(['Res. |fp| = ', num2str(norm(fp)/length(fp))]);
        fprintf( 'PH  iterations took: %2.4f s\n', toc )
        x = [u; p];
        min(x)
        max(x)

%         figure(11)
%         subplot(221)
%         spy(aJ)
%         subplot(222)
%         spy(bJ)
%         subplot(223)
%         spy(eJ)
%         subplot(224)
%         spy(fJ)
    end
end

eqp=reshape(eqp,nx-1,nz-1)';
eqv=reshape(eqv,nx+1,nz)';
equ=reshape(equ,nx,nz+1)';
eta_n=reshape(eta_n,nx-1,nz-1)';
eta_s=reshape(eta_s,nx,nz)';

tag_s=reshape(tag_s,nx,nz)';
tag_n=reshape(tag_n,nx-1,nz-1)';
tag_u=reshape(tag_u,nx,nz+1)';
tag_v=reshape(tag_v,nx+1,nz)';
 
% eqp(end-15:end, 1:7)
% eqv(end-15:end, 1:7)
% equ(end-15:end, 1:7)
% 
% eta_n(end-15:end, 1:7)
% eta_s(end-15:end, 1:7)
% 
% tag_s(end-15:end, 1:7)

% sum(isnan(x))
% find(isnan(x)==1)

% build a staggered grid with tags
if tags == 1
figure(2)
hold on
mesh(xg_coord, zg_coord, zeros(length(zg_coord),length(xg_coord)))

[xc2, zc2] = meshgrid(xc_coord, zc_coord);
xc2 = xc2(:);
zc2 = zc2(:);
tag_n = tag_n(:);
plot(xc2(tag_n==-1), zc2(tag_n==-1), 'or')
plot(xc2(tag_n==31), zc2(tag_n==31), 'og')


[xn2, zvx2] = meshgrid(xg_coord, zvx_coord);
xn2 = xn2(:);
zvx2 = zvx2(:);
tag_u = tag_u(:);
plot(xn2(tag_u==-1), zvx2(tag_u==-1), '+b', 'MarkerSize', 5)
plot(xn2(tag_u==11), zvx2(tag_u==11), '*b', 'MarkerSize', 5)
plot(xn2(tag_u==13), zvx2(tag_u==13), '*b', 'MarkerSize', 5)
plot(xn2(tag_u==0), zvx2(tag_u==0), '+b', 'MarkerSize', 5)

[xvz2, zn2] = meshgrid(xvz_coord, zg_coord);
xvz2 = xvz2(:);
zn2 = zn2(:);
tag_v = tag_v(:);
plot(xvz2(tag_v==-1), zn2(tag_v==-1), 'dk')
plot(xvz2(tag_v==11), zn2(tag_v==11), '*k')
plot(xvz2(tag_v==13), zn2(tag_v==13), '*k')
plot(xvz2(tag_v==0), zn2(tag_v==0), 'dk')

[xn2, zn2] = meshgrid(xg_coord, zg_coord);

plot(xn2(tag_s==-1), zn2(tag_s==-1), '.m')
% plot(xn2(tag_s==30), zn2(tag_s==30), '.y')
plot(xc2(tag_n==0), zc2(tag_n==0), 'dr')

if topo == 1
            plot( xtopo, ztopo, '.c' )
end
end

if visc == 1
    figure(22)
    subplot(211)
    imagesc(xg_coord, zg_coord, log10(eta_s));
    shading flat, axis xy image, colorbar;
    
    eta_n(isnan(eta_n)==1) = 1e10;
    subplot(212)
    imagesc(xc_coord, zc_coord, log10(eta_n));
    shading flat, axis xy image, colorbar;
end




% [a,b]= find(equ==9386)
% plot( xg_coord(b), zvx_coord(a), '*' )
% 
% [a,b]= find(equ==9387)
% plot( xg_coord(b), zvx_coord(a), '*' )


% axis([-2e4 0 -1e5 5e4])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% filename = ['~/Desktop/Reproducibility/STRESS_BC/DANI/r21/Matrix.gzip_4cpu.h5'];
% 
% Ic  = hdf5read(filename,'/matrix/I');
% J   = hdf5read(filename,'/matrix/J');
% Ic  = Ic+1;
% J   = J+1;
% V2  = hdf5read(filename,'/matrix/V');
% L   = hdf5read(filename,'/matrix/rhs');
% eta = hdf5read(filename,'/matrix/eta_cell');
% 
% equ = hdf5read(filename,'/numbering/eqn_u');
% eqv = hdf5read(filename,'/numbering/eqn_v');
% eqp = hdf5read(filename,'/numbering/eqn_p');
% 
% params = hdf5read(filename,'/model/params');
% 
% nx = params(1);
% nz = params(2);
% dx = params(3);
% dz = params(4);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Decompress index I
% I = zeros(size(J));
% 
% m = 0;
% 
% I(1) = 1;
% p    = 1;
% 
% for k=2:length(Ic)
%     
%     off = Ic(k) - Ic(k-1);
%     
%     for l=1:off
%         I(m+l) = p;
%     end
%     m = m + off;
%     p = p + 1;
%  
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% rhs2 = L;
% I2 = cast(I, 'double');
% J2 = cast(J, 'double');
% M2 = sparse(I2,J2,V2);
% max(J2)
% J2
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% 
% % subplot 121
% % spy(M)
% % 
% % subplot 122
% % spy(M2)
% 
% difference = abs(abs(M2-M));
% difference = difference(difference>1e13);
% difference = sparse(difference);
% spy(difference);


% sym = abs(M);
% sym(sym<1e-3) = 0;
% max(sym)
% spy(sym);
% rhs-rhs2



% %sol = M\L;
% sol = L*0;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % sol = M\L;
% e  = zeros(nz-1,nx-1);
% 
% min(min(M-(M-M')));
% 
% K = M(1:eqp(1),1:eqp(1));
% % K = M(eqp(1)+1:end,eqp(1)+1:end);
% KK=K-K';
% spy(KK)
% 
% % vx = zeros(nz+1,nx);
% % vz = zeros(nz,nx+1);
% % p  = zeros(nz-1,nx-1);
% % 
% % eu = reshape(equ,nx,nz+1)';
% % ev = reshape(eqv,nx+1,nz)';
% % ep = reshape(eqp,nx-1,nz-1)';
% % 
% % for j=1:nz+1
% %     for i =1:nx
% %         if ( (eu(j,i)) ~= -1)
% %             vx(j,i) = sol(eu(j,i)+1);
% %         end
% %     end
% % end
% % 
% % for j=1:nz
% %     for i =1:nx+1
% %         if (ev(j,i) ~= -1)
% %             vz(j,i) = sol(ev(j,i)+1);
% %         end
% %     end
% % end
% % 
% % for j=1:nz-1
% %     for i =1:nx-1
% %         if (ep(j,i) ~= -1)
% %             p(j,i) = sol(ep(j,i)+1);
% %         end
% %     end
% % end
% % p = p-mean(p(:));
% % M = full(M);
% % 
% % 
% % figure(2)
% % subplot(311)
% % imagesc(vx)
% % axis image
% % colorbar
% % 
% % subplot(312)
% % imagesc(vz)
% % axis image
% % colorbar
% % 
% % subplot(313)
% % imagesc(p)
% % axis image
% % colorbar
% 
% for j=1:nz-1
%     for i =1:nx-1
%         if (ep(j,i) ~= -1)
%             p(j,i) = sol(ep(j,i)+1);
%         end
%         cc = i + (j-1)*(nx-1); 
%         e(j,i) = eta(cc);
%     end
% end
% p = p-mean(p(:));
% 
% 
% 
% % Define Elman/Wathan style upper
% M_pc = [ K, G ; zeros(size(D)), -S_pc ];
% 
% % Define inverse of Elman/Wathan style upper
% %iiM_pc = [ iiK, G ; zeros(size(D)), -inv(S_pc) ];
% 
% % Right preconditioned operator
% %M_pc_operator = M * iiM_pc;
% % Left preconditioned operator
% % M_pc_operator = iiM_pc * M;
% % 
% % eeM = eig(full(M_pc_operator));
% % 
% % fprintf(1,'M: min(Re) %1.4e\n',min(real(eeM)));
% % fprintf(1,'M: max(Re) %1.4e\n',max(real(eeM)));
% % fprintf(1,'M: min(Im) %1.4e\n',min(imag(eeM)));
% % fprintf(1,'M: max(Im) %1.4e\n',max(imag(eeM)));
% % 
% % figure(6)
% % scatter(real(eeM),imag(eeM))
% 
% 
% 
% 
% %gmres(M,L,30,1.0e-5,100,M_pc);
% 
% UVP = gmres(M,L,30,1.0e-7,100,@(y)StokesPC_EWUpperTri(y,K,G,-S_pc,ev) );
% 
% rr = M*UVP - L;
% fprintf(1,'Norm2(rr) = %1.4e\n',norm(rr));
% fprintf(1,'Norm2(rr)/sq(N) = %1.4e\n',norm(rr)/sqrt(size(rr,1)));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % % plot 2D solution fields
% % figure(2)
% % subplot(311)
% % imagesc(vx)
% % axis image
% % colorbar
% % 
% % subplot(312)
% % imagesc(vz)
% % axis image
% % colorbar
% % 
% % subplot(313)
% % imagesc(p)
% % axis image
% % colorbar
% % 
% % % plot 2D viscosity fields
% % figure(3)
% % imagesc(log10(e))
% % axis image
% % colorbar
% % 
% % % Spy system blocks
% % figure(4)
% % subplot(221)
% % spy(K)
% % 
% % subplot(222)
% % spy(G)
% % 
% % subplot(223)
% % spy(D)
% % 
% % subplot(224)
% % spy(C)
% end
% 
% function y = StokesPC_EWUpperTri(r,K,G,S_pc,ev)
%     % split r => (ru,rp)
%     ru = r(1:max(ev(:)+1));   
%     rp = r(max(ev(:)+2):end);
%     
%     % apply factor
%     yp = S_pc\rp;    
%     yu = K\(ru - G*yp);
%        
%     % merge yu,yp => y
%     y = [yu;yp];
% end

