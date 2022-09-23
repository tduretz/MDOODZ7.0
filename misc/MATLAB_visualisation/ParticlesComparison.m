clear
istep = 0;
path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE/'

shift = 1;

filename_p0 = [path,'Particles',num2str(istep,'%05d'),'.gzip.h5'];
filename_p1 = [path,'Particles',num2str(istep,'%05d'),'_shift.gzip.h5'];

Part0.x     = hdf5read(filename_p0,'/Particles/x');
Part0.z     = hdf5read(filename_p0,'/Particles/z');
Part0.ph    = hdf5read(filename_p0,'/Particles/phase');
Part0.gen   = hdf5read(filename_p0,'/Particles/generation');
Part0.ind   = hdf5read(filename_p0,'/Particles/index');

Part1.x     = hdf5read(filename_p1,'/Particles/x');
Part1.z     = hdf5read(filename_p1,'/Particles/z');
Part1.ph    = hdf5read(filename_p1,'/Particles/phase');
Part1.gen   = hdf5read(filename_p1,'/Particles/generation');
Part1.ind   = hdf5read(filename_p1,'/Particles/index');

norm(Part0.x-Part1.x)
norm(Part0.z-(Part1.z+shift))