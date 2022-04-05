function main()
clear

% Compare initial pressure versus pressure at steady state

fsteady   = './CompFreeSurf/Output00020.gzip.h5';
[rsteady] = Read(fsteady);

fini   = './Output00001.gzip.h5';
[rini] = Read(fini);

Pini =  rini.P(:,fix(rini.nx/2))


dP = rsteady.P(:,fix(rsteady.nx/2)) - Pini;

figure(1), clf
subplot(131),hold on
plot( rini.Vz(:,fix(rini.nx/2)), rini.zg, '-b' )
plot( rsteady.Vz(:,fix(rini.nx/2)), rsteady.zg, '-r' )
subplot(132),hold on
% % imagesc(r51.xc, r51.zc, r51.P), colorbar, axis xy image
plot( Pini, rsteady.zc, '-k' )
plot( rini.P(:,fix(rini.nx/2)), rini.zc, '+b' )
title('Pressures')
subplot(133)
plot( dP, rini.zc, '-b' )
title('Pressure difference')
end

function [r51] = Read(fname51)
r51.p       = hdf5read(fname51, '/Model/Params');
r51.P       = hdf5read(fname51,'/Centers/P');
r51.P       = cast(r51.P , 'double');
r51.nx      = r51.p(4);
r51.nz      = r51.p(5);
r51.P       = reshape(r51.P,r51.nx-1,r51.nz-1)';
r51.z_grid  = hdf5read(fname51,'/Topo/z_grid');
r51.z_grid  = cast(r51.z_grid, 'double');
r51.xg      = hdf5read(fname51,'/Model/xg_coord');  
r51.xg      = cast(r51.xg, 'double');
r51.zg      = hdf5read(fname51,'/Model/zg_coord');  
r51.zg      = cast(r51.zg, 'double');
r51.xc      = hdf5read(fname51,'/Model/xc_coord');  
r51.xc      = cast(r51.xc, 'double');
r51.zc      = hdf5read(fname51,'/Model/zc_coord');  
r51.zc      = cast(r51.zc, 'double');
r51.Vx = hdf5read(fname51,'/VxNodes/Vx');
r51.Vz = hdf5read(fname51,'/VzNodes/Vz');
r51.Vx = cast(r51.Vx, 'double');
r51.Vz = cast(r51.Vz, 'double');
r51.Vx = reshape(r51.Vx,r51.nx,r51.nz+1)';
r51.Vz = reshape(r51.Vz,r51.nx+1,r51.nz)';
r51.dx = r51.xg(2) - r51.xg(1);
r51.dz = r51.zg(2) - r51.zg(1);
end