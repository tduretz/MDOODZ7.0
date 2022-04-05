% function Main

clear 
close all
clc

addpath '/Users/tduretz/SuiteSparse4.4.3/MATLAB_Tools/'
addpath '/Users/tduretz/SuiteSparse4.4.3/CHOLMOD/MATLAB/'
addpath '/Users/tduretz/SuiteSparse4.4.3/CXSparse/MATLAB/CSparse/'
addpath '/Users/tduretz/SuiteSparse4.4.3/UMFPACK/MATLAB/'
addpath '/Users/tduretz/SuiteSparse4.4.3/AMD/MATLAB/'
addpath '/Users/tduretz/SuiteSparse4.4.3/CAMD/MATLAB/'
addpath '/Users/tduretz/SuiteSparse4.4.3/CCOLAMD/MATLAB/'


% here you change this parameter to 16, 32 or 64 (resolution)
res = 21;

solve = 1;

istep = 1;

tags = 0;
visc = 0;
topo = 0;

jac = 0;

path = '/Users/tduretz/REPO_GIT/MDOODZ6.0/SOURCE/'


filename = [path 'Jacobian.gzip_1cpu.h5'];

IcA  = hdf5read(filename,'/matrix/IA');
JA   = hdf5read(filename,'/matrix/JA');
IcA  = IcA+1;
JA   = JA+1;
V1A  = hdf5read(filename,'/matrix/VA');

IcB  = hdf5read(filename,'/matrix/IB');
JB   = hdf5read(filename,'/matrix/JB');
IcB  = IcB+1;
JB   = JB+1;
V1B  = hdf5read(filename,'/matrix/VB');


IcC  = hdf5read(filename,'/matrix/IC');
JC   = hdf5read(filename,'/matrix/JC');
IcC  = IcC+1;
JC   = JC+1;
V1C  = hdf5read(filename,'/matrix/VC');

IcD  = hdf5read(filename,'/matrix/ID');
JD   = hdf5read(filename,'/matrix/JD');
IcD  = IcD+1;
JD   = JD+1;
V1D  = hdf5read(filename,'/matrix/VD');

rSt   = hdf5read(filename,'/matrix/rhs_mom');
rPt   = hdf5read(filename,'/matrix/rhs_cont');


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
IA = zeros(size(JA));

m = 0;

IA(1) = 1;
p    = 1;

for k=2:length(IcA)
    
    off = IcA(k) - IcA(k-1);
    
    for l=1:off
        IA(m+l) = p;
    end
    m = m + off;
    p = p + 1;
 
end

I1A  = cast(IA, 'double');
J1A  = cast(JA, 'double');
MA    = sparse(I1A, J1A, V1A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompress index I
IB = zeros(size(JB));

m = 0;

IB(1) = 1;
p    = 1;

for k=2:length(IcB)
    
    off = IcB(k) - IcB(k-1);
    
    for l=1:off
        IB(m+l) = p;
    end
    m = m + off;
    p = p + 1;
 
end

I1B  = cast(IB, 'double');
J1B  = cast(JB, 'double');
MB    = sparse(I1B, J1B, V1B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompress index I
IC = zeros(size(JC));

m = 0;

IC(1) = 1;
p    = 1;

for k=2:length(IcC)
    
    off = IcC(k) - IcC(k-1);
    
    for l=1:off
        IC(m+l) = p;
    end
    m = m + off;
    p = p + 1;
 
end

I1C  = cast(IC, 'double');
J1C  = cast(JC, 'double');
MC    = sparse(I1C, J1C, V1C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decompress index I
ID = zeros(size(JD));

m = 0;

ID(1) = 1;
p    = 1;

for k=2:length(IcD)
    
    off = IcD(k) - IcD(k-1);
    
    for l=1:off
        ID(m+l) = p;
    end
    m = m + off;
    p = p + 1;
 
end

I1D  = cast(ID, 'double');
J1D  = cast(JD, 'double');
MD    = sparse(I1D, J1D, V1D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
subplot 221
spy(MA)
subplot 222
spy(MB)
subplot 223
spy(MC)
subplot 224
spy(MD)

save('Jacobian','MA', 'MB', 'MC', 'MD')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2),
imagesc(reshape(eta_s,nx,nz)), colorbar

figure(3)
muc  = 3.2095e+22
filename = '/Users/tduretz/REPO_GIT/MDOODZ6.0/SOURCE/Output00000.gzip.h5'


       eta_n = hdf5read(filename,'/Centers/eta_n');  eta_n = cast(eta_n, 'double');
        eta_n   = reshape(eta_n,nx-1,nz-1)';
        
        
       eta_s = hdf5read(filename,'/Vertices/eta_s');  eta_s = cast(eta_s, 'double');
        eta_s   = reshape(eta_s,nx,nz)';

    imagesc((eta_s/muc)'), colorbar, set (gca,'ydir','normal')
