
function main
clear

path = './';

istep = 50;

file1 = [path,'Output',num2str(istep,'%05d'),'_ref.gzip.h5'];
file2 = [path,'Output',num2str(istep,'%05d'),'_low.gzip.h5'];
file3 = [path,'Output',num2str(istep,'%05d'),'_up.gzip.h5'];


[f1] = Read(file1);
[f2] = Read(file2);
[f3] = Read(file3);

figure(1), clf
% Viscosity
% subplot(131), imagesc(f1.xc_coord,f1.zc_coord,f1.eta_n), axis image xy, colorbar, shading flat
% subplot(132), imagesc(f1.xc_coord,f1.zc_coord,f1.eta_n), axis image xy, colorbar, shading flat
% subplot(133), imagesc(f3.xc_coord,f3.zc_coord,f3.eta_n), axis image xy, colorbar, shading flat

subplot(131), imagesc(f1.xc_coord,f1.zc_coord,f1.P), axis image xy, colorbar, shading flat
subplot(132), imagesc(f1.xc_coord,f1.zc_coord,f1.P), axis image xy, colorbar, shading flat
subplot(133), imagesc(f3.xc_coord,f3.zc_coord,f3.P), axis image xy, colorbar, shading flat

disp('viscosity')
norm(f1.eta_s-f2.eta_s)
norm(f1.eta_s-f3.eta_s)
norm(f1.eta_n-f2.eta_n)
norm(f1.eta_n-f3.eta_n)

disp('pressure')
norm(f1.P - f2.P)
norm(f1.P - f3.P)

disp('velocity')
norm(f1.Vx - f2.Vx)
norm(f1.Vx - f3.Vx)
norm(f1.Vz - f2.Vz)
norm(f1.Vz - f3.Vz)

% ----------------------------------------- %
end

function [f1] = Read(file1)

f1.params     = hdf5read(file1,'/Model/Params');
f1.xg_coord   = hdf5read(file1,'/Model/xg_coord');  f1.xg_coord  = cast(f1.xg_coord, 'double');
f1.zg_coord   = hdf5read(file1,'/Model/zg_coord');  f1.zg_coord  = cast(f1.zg_coord, 'double');
f1.xc_coord   = hdf5read(file1,'/Model/xc_coord');  f1.xc_coord  = cast(f1.xc_coord, 'double');
f1.zc_coord   = hdf5read(file1,'/Model/zc_coord');  f1.zc_coord  = cast(f1.zc_coord, 'double');
f1.xvz_coord  = hdf5read(file1,'/Model/xvz_coord'); f1.xvz_coord = cast(f1.xvz_coord, 'double');
f1.zvx_coord  = hdf5read(file1,'/Model/zvx_coord'); f1.zvx_coord = cast(f1.zvx_coord, 'double');
f1.VizGrid.x  = hdf5read(file1,'/VizGrid/xviz');    f1.VizGrid.x = cast(f1.VizGrid.x, 'double');
f1.VizGrid.z  = hdf5read(file1,'/VizGrid/zviz');    f1.VizGrid.z = cast(f1.VizGrid.z, 'double');
f1.eta_n = hdf5read(file1,'/Centers/eta_n');  f1.eta_n = cast(f1.eta_n, 'double');
f1.eta_n = reshape(f1.eta_n,f1.params(4)-1,f1.params(5)-1)';
f1.eta_s = hdf5read(file1,'/Vertices/eta_s');  f1.eta_s = cast(f1.eta_s, 'double');
f1.eta_s = reshape(f1.eta_s,f1.params(4)-0,f1.params(5)-0)';
f1.P     = hdf5read(file1,'/Centers/P'); f1.P  = cast(f1.P , 'double'); 
f1.P     = reshape(f1.P,f1.params(4)-1,f1.params(5)-1)';
f1.Vx    = hdf5read(file1,'/VxNodes/Vx'); f1.Vx = cast(f1.Vx, 'double');
f1.Vz    = hdf5read(file1,'/VzNodes/Vz'); f1.Vz = cast(f1.Vz, 'double');
f1.Vx    = reshape(f1.Vx,f1.params(4)+0,f1.params(5)+1)';
f1.Vz    = reshape(f1.Vz,f1.params(4)+1,f1.params(5)+0)';

end

