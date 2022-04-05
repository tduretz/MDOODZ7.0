clear

fname51   = './Clement_51x161/Output00050.gzip.h5';
fname101  = './Clement_101x161/Output00050.gzip.h5';
fname201  = './Clement_201x161/Output00050.gzip.h5';
fname501  = './Clement_501x161/Output00050.gzip.h5';
fname1001 = './Clement_1001x161/Output00082.gzip.h5';
fname2001 = './Clement_2001x161/Output00150.gzip.h5';

r51.p   = hdf5read(fname51, '/Model/Params');
r101.p  = hdf5read(fname101, '/Model/Params');
r201.p  = hdf5read(fname201, '/Model/Params');
r501.p  = hdf5read(fname501, '/Model/Params');
r1001.p = hdf5read(fname1001,'/Model/Params');
r2001.p = hdf5read(fname2001,'/Model/Params');

r51.z_grid  = hdf5read(fname51,'/Topo/z_grid');
r51.z_grid  = cast(r51.z_grid, 'double');
r51.xg      = hdf5read(fname51,'/Model/xg_coord');  
r51.xg      = cast(r51.xg, 'double');

r101.z_grid  = hdf5read(fname101,'/Topo/z_grid');
r101.z_grid  = cast(r101.z_grid, 'double');
r101.xg      = hdf5read(fname101,'/Model/xg_coord');  
r101.xg      = cast(r101.xg, 'double');

r201.z_grid  = hdf5read(fname201,'/Topo/z_grid');
r201.z_grid  = cast(r201.z_grid, 'double');
r201.xg      = hdf5read(fname201,'/Model/xg_coord');  
r201.xg      = cast(r201.xg, 'double');

r501.z_grid  = hdf5read(fname501,'/Topo/z_grid');
r501.z_grid  = cast(r501.z_grid, 'double');
r501.xg      = hdf5read(fname501,'/Model/xg_coord');  
r501.xg      = cast(r501.xg, 'double');

r1001.z_grid  = hdf5read(fname1001,'/Topo/z_grid');
r1001.z_grid  = cast(r1001.z_grid, 'double');
r1001.xg      = hdf5read(fname1001,'/Model/xg_coord');  
r1001.xg      = cast(r1001.xg, 'double');

r2001.z_grid  = hdf5read(fname2001,'/Topo/z_grid');
r2001.z_grid  = cast(r2001.z_grid, 'double');
r2001.xg      = hdf5read(fname2001,'/Model/xg_coord');  
r2001.xg      = cast(r2001.xg, 'double');

r51.p(1)
r101.p(1)
r201.p(1)
r501.p(1)
r1001.p(1)
r2001.p(1)

figure(1), clf
hold on
plot(r51.xg, r51.z_grid, '-k')
plot(r101.xg, r101.z_grid, '.-k')
plot(r201.xg, r201.z_grid, '-.k')
plot(r501.xg, r501.z_grid, 'g-')
plot(r1001.xg, r1001.z_grid, '-b')
plot(r2001.xg, r2001.z_grid, '-r')