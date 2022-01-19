function [ MA, MB, MC, MD, rSt, rPt, fSt, fPt ] = ExtractDecoupledMatrixFromMDoodz( filename )

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

fSt   = hdf5read(filename,'/matrix/F_mom');
fPt   = hdf5read(filename,'/matrix/F_cont');

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

end

