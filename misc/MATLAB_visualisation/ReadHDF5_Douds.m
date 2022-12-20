
function Main

clear
clear all
% close all
clc
<<<<<<< HEAD
path = '~/REPO_GIT/MDOODZ7.0/MDLIB/';
=======

DEBUG = 0;

MarkSize=1e0 ;

mdoodz6       = 1;

path = '/Users/tduretz/REPO_GIT/MDOODZ6.0/SOURCE/Compression/';
path = '/Users/tduretz/REPO_GIT/MDOODZ6.0/SOURCE/Shear_periodic_VEVP/';
path = '/Users/tduretz/REPO_GIT/MDOODZ6.0/SOURCE/';
path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE/Huismans_Reg2e20_noise/';
path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE/Huismans_OS2e5_n1_5/';
path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE/Huismans_Reg8e19_noise_NEW/';

path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE//';


% path = '/Users/imac/REPO_GIT/mdoodz4.5_repo/SOURCE/'
% mdoodz6       = 0;


% path = '/Volumes/Seagate4TB/Wedge_MD6/LR/'
% path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE/RUN_Ridge_40/'
% % path = '/Volumes/Seagate4TB/LithoScale/Ext_HR2_5e20/'
% % path = '/Users/imac/REPO_GIT/MDOODZ6.0_beauty/SOURCE/'
% path = '/Volumes/Seagate4TB/LithoScale/N3/Comp_HHR_2e20_Courant0.25/'
% path = '/Volumes/Seagate4TB/LithoScale/N3/Comp_HHR_2e20_Courant0.25/'
% path = '/Volumes/Seagate4TB/LithoScale/N3/Comp_HR_4e20/'

% path = '/Users/imac/Downloads/MDoodz4.5_Plateau/'
% mdoodz6       = 0;

% path = '/Volumes/Passport/test_AL/G_1e12/'

>>>>>>> main
cd(path)
DEBUG       = 0;
mdoodz6     = 1;

% Files
istart = 10;
ijump  = 1;
iend   = 10;

%--------------------------------------------------
% what do you want to plot:
%--------------------------------------------------
eta_sym         = 0;
eta_plot        = 0;
rho_plot        = 0;
phase_on_grid   = 0;
phase_temp2     = 0;
vel_plot        = 0;
vel_vectors     = 0;
vel_divergence  = 0;
pre_plot        = 0;
dyna_pre        = 0;
stress_inv      = 1;
stress_evol     = 0;
stress_plot     = 0;
srate_plot      = 0;
acc_strain      = 0;
temperature     = 0;
temp_evol       = 0;
Stefan          = 0;
paper           = 0;
kinetic_e       = 0;
strainR_div     = 0;
paper2          = 0;
Christmas_Tree  = 0;
topo            = 0;
topo_eta_plot   = 0;
topo_SR_plot    = 0;
topo_maps       = 0;
phases_uplift   = 0;
dt_time         = 0;
srate_add       = 1;
acc_strain_add  = 0;
phase_temp      = 0;
srate_add_perc  = 0;
grain_size      = 0;
grain_size_evol = 0;
gs_stefan       = 0;
X_plot          = 0;
srate_ratio     = 0;
PTmax           = 0;
TminusT0        = 0;
emergency_benoit = 0;
fstrain          = 0;
shear_heating    = 0;
princi_stress    = 0;
director_vector  = 0;
Pl_soft          = 0;
Sole             = 0;
overstress       = 0;

% Visualisation options
printfig      = 0;
print2screen  = 1;
res           = '-r300';
format        = '-dpng';
color_map     = 'Benoit_MLPS';
color_map     = 'geoffroy_ast3';
color_map     = 'JPB2';
Ftsz          = 20;
file_suffix   = '';
ConvTest      = 0;
show          = 0;
Ccontours     = 0;
step          = 10;
MaskAir       = 0;
ColorFabio    = 1;
y             = 365.25*24*3600;
My            = 1e6*y;
topo_ref      = 0; %6370

% minVx = -2.5e-13;
% maxVx =  2.5e-13;
% minVz = -1.5e-13;
% maxVz =  0.0e-13;

% Contour
eps_val = 5e-14;
% eps_val = 1e4;

% minRho = 3000;
% maxRho = 3600;

minX = 0;
maxX = 1;

% Colorbar limits








% minHkm = -0.15; 
% maxHkm = 0.1; % for comparison

% minVxmmy = -0.1;
% maxVxmmy = 0.1; % for comparison
% 
% minVymmy = 0;
% maxVymmy = 1; % for comparison


% minVxmmy = -10;
% maxVxmmy =  10;  % for PaulineRift
% minVymmy = -2.5;
% maxVymmy =  2.5;  % for PaulineRift


% % minHkm = -1.5; 
% % maxHkm = 0.5; % for comparison
% % minVxmmy = -2.5;
% % maxVxmmy =  2.5; % for comparison
% % minVymmy = -5;
% % maxVymmy = 15; % for comparison

minTxx = -11e0;
maxTxx =  0.642;

minStr = 0.01;
maxStr = 0.4;

minPdyn = -1e8;
maxPdyn =  5e8;

% minEta = 19;
% maxEta = 25;

mindiv  =-0.25e-14;
maxdiv  = 0.25e-14;

% 
% minEii = -17;
% maxEii = -15;
% 
% minSii = 1e6;
% maxSii = 9e7;

% minEii = -17;
% maxEii = -13;
% 
% minSii = 1e6;
% maxSii = 400e6;

% Size of the window
crop       = 0;

lim.xmin   = -4;
lim.xmax   =  4;
lim.zmin   = -4;
lim.zmax   =  1;

% lim.xmin   = -75;
% lim.xmax   =  75;
% lim.zmin   = -35;
% lim.zmax   = 5;

% lim.xmin   = 45;
% lim.xmax   =  70;
% lim.zmin   = -6.0;
% lim.zmax   = 2;

% minP = 0;
% maxP = 1e9;

% In order to calculate lithostatic pressure
gz = -9.81;

% Arrays
dt_vec = zeros(10000,1);
E_kin  = zeros(10000,1);
t_vec  = zeros(10000,1);
gamma  = zeros(10000,1);
error  = zeros(10000,1);
etapp  = zeros(10000,1);
mxtopo = zeros(10000,1);
icount = 0;
time   = 0;
Ws     = 0;
Es     = 0;

SbThickYuri = zeros( 10000,1);
SbThick     = zeros( 10000,1);
SbDT        = zeros( 10000,1);
SbDTY       = zeros( 10000,1);
SbAngle     = zeros( 10000,1);
SbEiiFact   = zeros( 10000,1);
short       = zeros( 10000,1);
siivec      = zeros( 10000,1);
eiivec      = zeros( 10000,1);
eiimaxvec   = zeros( 10000,1);
T_postrift  = zeros( 10000,1);
T_prerift   = zeros( 10000,1);
timevec     = zeros( 10000,1);

strainvec    = zeros( 10000,1);
neck_thick   = zeros( 10000,1);
swell_thick  = zeros( 10000,1);
neck_stress  = zeros( 10000,1);
swell_stress = zeros( 10000,1);
neck_gs      = zeros( 10000,1);
swell_gs     = zeros( 10000,1);
neck_srate   = zeros( 10000,1);
swell_srate  = zeros( 10000,1);
neck_T       = zeros( 10000,1);
swell_T      = zeros( 10000,1);

%--------------------------------------------------
% Define RGB colormaps for composition figures:
%--------------------------------------------------

if strcmp(color_map, 'inc2boules') == 1
    %             hold on, plot(Part.x(Part.ph==0),Part.z(Part.ph==0),'.','MarkerSize',MarkSize,'color',[0.028 0.54 0.295]);
    %             hold on, plot(Part.x(Part.ph==3),Part.z(Part.ph==3),'k.','MarkerSize',MarkSize);
    %             hold on, plot(Part.x(Part.ph==1),Part.z(Part.ph==1),'r.','MarkerSize',MarkSize);
    %             hold on, plot(Part.x(Part.ph==2),Part.z(Part.ph==2),'w.','MarkerSize',MarkSize);
elseif strcmp(color_map, 'piston') == 1
    map = [56/255 117/255 215/255; ...
        1/255 31/255 113/255;...
        207/255 21/255 16/255;...
        248/255 216/255 71/255;...
        237/255 237/255 237/255;...
        1 1 1;...
        237/255 237/255 237/255;....
        0 0 0;];
elseif strcmp(color_map, 'geoffroy') == 1
    map = [        1           1           1
        215/255      46/255     39/255;...
        250/255     159/255    130/255;...
        163/255     228/255    171/255;...
        29/255     168/255    158/255];
    
elseif strcmp(color_map, 'geoffroy_ast') == 1
    map = [        1           1           1
        215/255      46/255     39/255;...
        250/255     159/255    130/255;...
        163/255     228/255    171/255;...
        29/255     168/255    158/255;...
        90/255     190/255    160/255];
    
elseif strcmp(color_map, 'geoffroy_ast2') == 1
    map = [        255           255           255
        232 53 54;
        255 242 164;...
        240 153 126;...
        248 207 189;...
        235 51 46;...
        201 143 122;...
        
        143 198 153;...
        101 139 110 ;...
        213 229 208]./255;
    
elseif strcmp(color_map, 'geoffroy_ast3') == 1
    map = [        255           255           255
        232 53 54;
        247 207 190;
        165 227 173;
        98 117 103;
        213 228 209;
        200 143 124]./255;
elseif strcmp(color_map, 'JPB') == 1
    map = [        255           255           255
        105 194 239;
        247 207 190;
        221 113  43;
        232 53 54;
        196 140 124;
        57 145  47;
        183 222 133;
        155 154 229;
        055 054 129;]./255;
elseif strcmp(color_map, 'JPB2') == 1
    map = [        255           255           255
        105 194 239;
        247 207 190;
        221 113  43;
        232 53 54;
        196 140 124;
        57 145  47;
        183 222 133;
        155 154 229;
        055 054 129;
        239, 139, 176;]./255;
    
elseif strcmp(color_map, 'Benoit_MLPS') == 1
    map = 'jet';
else
    map = 'jet';
    
end

% Create directories
if printfig == 1 || printfig == 2
    
    if (exist('./Fig_Viscosity', 'dir') == 0)
        mkdir ./Fig_Viscosity
    end
    
    if (exist('./Fig_Density', 'dir') == 0)
        mkdir ./Fig_Density
    end
    
    if (exist('./Fig_Phases', 'dir') == 0)
        mkdir ./Fig_Phases
    end
    
    if (exist('./Fig_Velocity', 'dir') == 0)
        mkdir ./Fig_Velocity
    end
    
    if (exist('./Fig_Pressure', 'dir') == 0)
        mkdir ./Fig_Pressure
    end
    
    if (exist('./Fig_Stress_StrainRate', 'dir') == 0)
        mkdir ./Fig_Stress_StrainRate
    end
    
    if (exist('./Fig_StressComponents', 'dir') == 0)
        mkdir ./Fig_StressComponents
    end
    
    if (exist('./Fig_StrainRateComponents', 'dir') == 0)
        mkdir ./Fig_StrainRateComponents
    end
    
    if (exist('./Fig_VelocityMagnitude', 'dir') == 0)
        mkdir ./Fig_VelocityMagnitude
    end
    
    if (exist('./Fig_VelocityDivergence', 'dir') == 0)
        mkdir ./Fig_VelocityDivergence
    end
    
    if (exist('./Fig_DynamicPressure', 'dir') == 0)
        mkdir ./Fig_DynamicPressure
    end
    
    if (exist('./Fig_AccStrain', 'dir') == 0)
        mkdir ./Fig_AccStrain
    end
    
    if (exist('./Fig_Temperature', 'dir') == 0)
        mkdir ./Fig_Temperature
    end
    
    if (exist('./Fig_GS', 'dir') == 0)
        mkdir ./Fig_GS
    end
    
    if (exist('./Fig_Phi', 'dir') == 0)
        mkdir ./Fig_Phi
    end
    
    if (exist('./Fig_X', 'dir') == 0)
        mkdir ./Fig_X
    end
    
    if (exist('./Fig_Topo', 'dir') == 0)
        mkdir ./Fig_Topo
    end
    
    if (exist('./Fig_OverS', 'dir') == 0)
        mkdir ./Fig_OverS
    end
end

%% LOOP ON OUTPUT FILES
for istep=istart:ijump:iend
    
    istep
    
    tic
    
    figCount = 0;
    
    if DEBUG == 0
        filename = [path,'Output',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];
    else
        filename = [path,'DEBUG' ,num2str(istep,'%05d'),file_suffix,'.gzip.h5'];
    end
    
    if exist(filename,'file')~=0
        icount = icount + 1;
        
        params     = hdf5read(filename,'/Model/Params');
        xg_coord   = hdf5read(filename,'/Model/xg_coord');  xg_coord  = cast(xg_coord, 'double');
        zg_coord   = hdf5read(filename,'/Model/zg_coord');  zg_coord  = cast(zg_coord, 'double');
        xc_coord   = hdf5read(filename,'/Model/xc_coord');  xc_coord  = cast(xc_coord, 'double');
        zc_coord   = hdf5read(filename,'/Model/zc_coord');  zc_coord  = cast(zc_coord, 'double');
        xvz_coord  = hdf5read(filename,'/Model/xvz_coord'); xvz_coord = cast(xvz_coord, 'double');
        zvx_coord  = hdf5read(filename,'/Model/zvx_coord'); zvx_coord = cast(zvx_coord, 'double');
        VizGrid.x  = hdf5read(filename,'/VizGrid/xviz');    VizGrid.x = cast(VizGrid.x, 'double');
        VizGrid.z  = hdf5read(filename,'/VizGrid/zviz');    VizGrid.z = cast(VizGrid.z, 'double');

        try
            VizGrid.x_hr  = hdf5read(filename,'/VizGrid/xviz_hr');    VizGrid.x_hr = cast(VizGrid.x_hr, 'double');
            VizGrid.z_hr  = hdf5read(filename,'/VizGrid/zviz_hr');    VizGrid.z_hr = cast(VizGrid.z_hr, 'double');
        catch
            VizGrid.x_hr  = 0;
            VizGrid.z_hr  = 0;
        end
        if length(VizGrid.x) == length(xg_coord) || length(VizGrid.x) == (2*(length(xg_coord)-1)+1)
        VizGrid.x  = 0.5*(VizGrid.x(2:end)+VizGrid.x(1:end-1));
        VizGrid.z  = 0.5*(VizGrid.z(2:end)+VizGrid.z(1:end-1));
        try
        VizGrid.x_hr  = 0.5*(VizGrid.x_hr(2:end)+VizGrid.x_hr(1:end-1));
        VizGrid.z_hr  = 0.5*(VizGrid.z_hr(2:end)+VizGrid.z_hr(1:end-1));
        catch
            VizGrid.x_hr  = 0;
            VizGrid.z_hr  = 0;
        end
    end
        
        nx    = params(4);
        nz    = params(5);
        ncx   = nx-1;
        ncz   = nz-1;
        time0 = time;
        time  = params(1);
        t_vec(icount) = time;
        
        if params(1) < 60
            TimeLabel = [' t = ' num2str(time) ' s '];
        elseif time >= 60 && time < 3600
            TimeLabel = [' t = ' num2str(time/60) ' min '];
        elseif time >= 3600 && time < 24*3600
            TimeLabel = [' t = ' num2str(time/3600) ' h '];
        elseif time >= 24*3600 && time < 365.25*24*3600
            TimeLabel = [' t = ' num2str(params(1)/24/3600) ' j '];
        elseif time >= 365.25*24*3600 && time < 1e3*365.25*24*3600
            TimeLabel = [' t = ' num2str(params(1)/365/24/3600) ' a '];
        elseif time >= 1e3*365.25*24*3600 && time < 1e6*365.25*24*3600
            TimeLabel = [' t = ' num2str(time/1e3/365.25/24/3600) ' ka '];
        elseif time >= 1e6*365.25*24*3600
            TimeLabel = [' t = ' num2str(time/1e6/365.25/24/3600) ' Ma '];
        end
        
        
        % Spacing
        L = xg_coord(end) - xg_coord(1);
        H = zg_coord(end) - zg_coord(1);
        vol = H*L;
        dx = xg_coord(2) - xg_coord(1);
        dz = zg_coord(2) - zg_coord(1);
        
        if L < 1e-6
            
            xg_plot   = xg_coord/1e-0;
            zg_plot   = zg_coord/1e-0;
            xc_plot   = xc_coord/1e-0;
            zc_plot   = zc_coord/1e-0;
            xvz_plot  = xvz_coord/1e-0;
            zvx_plot  = zvx_coord/1e-0;
            VizGrid.x_plot = VizGrid.x/1e-0;
            VizGrid.z_plot = VizGrid.z/1e-0;
            VizGrid.x_plot_hr = VizGrid.x_hr/1e-0;
            VizGrid.z_plot_hr = VizGrid.z_hr/1e-0;
            xLabel = ['x [nm]'];
            zLabel = ['y [nm]'];
            
        elseif L < 1e3
            xg_plot   = xg_coord/1e0;
            zg_plot   = zg_coord/1e0;
            xc_plot   = xc_coord/1e0;
            zc_plot   = zc_coord/1e0;
            xvz_plot  = xvz_coord/1e0;
            zvx_plot  = zvx_coord/1e0;
            VizGrid.x_plot = VizGrid.x/1e0;
            VizGrid.z_plot = VizGrid.z/1e0;
            VizGrid.x_plot_hr = VizGrid.x_hr/1e-0;
            VizGrid.z_plot_hr = VizGrid.z_hr/1e-0;
            xLabel = ['x [m]'];
            zLabel = ['y [m]'];
        else
            xg_plot   = xg_coord/1e3;
            zg_plot   = zg_coord/1e3;
            xc_plot   = xc_coord/1e3;
            zc_plot   = zc_coord/1e3;
            xvz_plot  = xvz_coord/1e3;
            zvx_plot  = zvx_coord/1e3;
            VizGrid.x_plot = VizGrid.x/1e3;
            VizGrid.z_plot = VizGrid.z/1e3;
            VizGrid.x_plot_hr = VizGrid.x_hr/1e3;
            VizGrid.z_plot_hr = VizGrid.z_hr/1e3;
            xLabel = ['x [km]'];
            zLabel = ['y [km]'];
        end
        
        % Adaptive min max axis
        dt =  time - time0; 
        
        % Initial model length
        if istep == 0 && icount == 1
            L0 = xg_coord(end)-xg_coord(1);
            L_t = 0;
            lim0.xmin = xg_coord(1);
            lim0.xmax = xg_coord(end);
            lim0.zmin = 0;
            lim0.zmax = 1.5*lim0.xmax;
        elseif istep > 0 && icount == 1
            file0 = [path,'Output',num2str(0,'%05d'),file_suffix,'.gzip.h5'];
            if exist( file0, 'file') > 0
                xg_0  = hdf5read(file0,'/Model/xg_coord');
                L0 = xg_0(end)-xg_0(1);
            else
                L0 = xg_coord(end)-xg_coord(1);
            end
        end
        
        H = zg_coord(end)-zg_coord(1);
        strain    = (L0-L)/L0*100;
        
        if istep > 0 && icount>1
            %         OldFile = [path,'Output',num2str(istep-ijump,'%05d'),file_suffix,'.gzip.h5'];
            VX = hdf5read(filename,'/VxNodes/Vx'); VX = cast(VX, 'double');
            VX = reshape(VX,nx,nz+1)';
            VXBC = abs(mean(VX(:,1)));
            L_t = L_t + 2.0*VXBC*dt;
            strain2 = 100-abs(L0-L_t)/L0*100;
            strain    = abs(L0-L)/L0*100;
            
        elseif icount==1
            L_t = 0;
            strain2 = 0;
            strain    = abs(L0-L)/L0*100;
            strain2 = strain;
        end
        
        if strain2 > 1e-13
            TimeLabel = [TimeLabel, ' - strain = ', num2str(strain), ' %'];
        end
        
%         % Generate a mask for inactive air cells over the free surface
%         eta = hdf5read(filename,'/Centers/eta_n');  eta = cast(eta, 'double');
%         eta   = reshape(eta,params(4)-1,params(5)-1)';
%         eta(eta==0) = nan;
%         [xc2,zc2]   = meshgrid(xc_plot, zc_plot);
%         x_air       =  xc2(isnan(eta));
%         z_air       =  zc2(isnan(eta));
%         x_tab       = zeros(length(x_air), 4);
%         z_tab       = zeros(length(x_air), 4);
%         f           = zeros(length(x_air),1);
%         id          = [ 1 2 4 3 ];
%         x_arr = [x_air(:) x_air(:) x_air(:)+dx/1e3 x_air(:)+dx/1e3];
%         z_arr = [z_air(:) z_air(:)+dz/1e3 z_air(:) z_air(:)+dz/1e3];
%         x_tab(:,:) = x_arr;
%         z_tab(:,:) = z_arr;
        
        % Generate a mask for inactive air cells over the free surface
        eta = hdf5read(filename,'/Vertices/eta_s');  eta = cast(eta, 'double');
        eta   = reshape(eta,params(4)-0,params(5)-0)';
        eta(eta==0) = nan;
        [xg2,zg2]   = meshgrid(xg_plot, zg_plot);
        x_air       =  xg2(isnan(eta));
        z_air       =  zg2(isnan(eta));
        x_tab       = zeros(length(x_air), 4);
        z_tab       = zeros(length(x_air), 4);
        f           = zeros(length(x_air),1);
        id          = [ 1 2 4 3 ];
        dxa         = dx/2/1e3;
        dza         = dz/2/1e3;
        x_arr = [x_air(:)-dxa x_air(:)-dxa x_air(:)+dxa x_air(:)+dxa];
        z_arr = [z_air(:)-dza z_air(:)+dza z_air(:)-dza z_air(:)+dza];
        x_tab(:,:) = x_arr;
        z_tab(:,:) = z_arr;
        
        %--------------------------------------------------
        % viscosity symmetry plot
        %--------------------------------------------------
        if (eta_sym==1)
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta_n = hdf5read(filename,'/Centers/eta_n');  eta_n = cast(eta_n, 'double');
            eta_s  = reshape(eta_s,params(4),params(5))';
            eta_n = reshape(eta_n,params(4)-1,params(5)-1)';
            nx         = size(eta_s,2);
            
            if mod(nx,2)==0
                eta_sym_s  = eta_s;
                eta_sym_s(:,1:nx/2    )  = abs(eta_sym_s(:,1:nx/2    ) -  eta_s(:,end:-1:nx/2+1));
                eta_sym_s(:,nx/2+1:end)  = abs(eta_sym_s(:,nx/2+1:end) -  eta_s(:,nx/2:-1:1    ));
                eta_sym_n  = eta_n;
                eta_sym_n(:,1:nx/2    )  = abs(eta_sym_n(:,1:nx/2    ) -  eta_n(:,end:-1:nx/2));
                eta_sym_n(:,nx/2+1:end)  = abs(eta_sym_n(:,nx/2+1:end) -  eta_n(:,nx/2-1:-1:1    ));
            else
                eta_sym_n  = eta_n;
                nx1 = nx-1;
                eta_sym_n(:,1:nx1/2    )  = abs(eta_sym_n(:,1:nx1/2    ) -  eta_n(:,end:-1:nx1/2+1));
                eta_sym_n(:,nx1/2+1:end)  = abs(eta_sym_n(:,nx1/2+1:end) -  eta_n(:,nx1/2:-1:1    ));
                eta_sym_s  = eta_s;
                nx1 = nx+1;
                eta_sym_s(:,1:nx1/2    )  = abs(eta_sym_s(:,1:nx1/2    ) -  eta_s(:,end:-1:nx1/2));
                eta_sym_s(:,nx1/2+1:end)  = abs(eta_sym_s(:,nx1/2+1:end) -  eta_s(:,nx1/2-1:-1:1    ));
            end
            
            % Nodes
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount); clf;
            else
                figure('Visible', 'Off');
            end
            
            subplot(221)
            imagesc(xg_plot, zg_plot, log10(eta_sym_s));
            shading flat, axis xy image, colorbar;
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            hold off
            %         title(['Viscosity (log10) at' TimeLabel])
            title(['eta s min = ', num2str(min(eta_sym_s(eta_sym_s>0)), '%2.4e'), ' max = ', num2str(max(eta_sym_s(eta_sym_s>0)), '%2.4e')])
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1, xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(223)
            imagesc(xc_plot, zc_plot, log10(eta_sym_n));
            shading flat, axis xy image, colorbar;
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            hold off
            %         title(['Viscosity (log10) at' TimeLabel])
            title(['eta n min = ', num2str(min(eta_sym_n(eta_sym_n>0)), '%2.4e'), ' max = ', num2str(max(eta_sym_n(eta_sym_n>0)), '%2.4e')])
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(222)
            imagesc(xg_plot, zg_plot, log10(eta_s));
            shading flat, axis xy image, colorbar;
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            hold off
            %         title(['Viscosity (log10) at' TimeLabel])
            title(['eta s min = ', num2str(min(eta_s(eta_s>0)), '%2.4e'), ' max = ', num2str(max(eta_s(eta_s>0)), '%2.4e')])
            xlabel(xLabel), ylabel(zLabel);
            if exist('minEta', 'var') caxis([minEta maxEta]); end
            if crop == 1, xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(224)
            imagesc(xc_plot, zc_plot, log10(eta_n));
            shading flat, axis xy image, colorbar;
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            hold off
            %         title(['Viscosity (log10) at' TimeLabel])
            title(['eta n min = ', num2str(min(eta_n(eta_n>0)), '%2.4e'), ' max = ', num2str(max(eta_n(eta_n>0)), '%2.4e')])
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minEta', 'var') caxis([minEta maxEta]); end
            drawnow
 
 
        end
        
        %--------------------------------------------------
        % viscosity plot
        %--------------------------------------------------
        if (eta_plot==1)
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta_n = hdf5read(filename,'/Centers/eta_n');  eta_n = cast(eta_n, 'double');
            eta   = reshape(eta_s,params(4),params(5))';
            eta_n = reshape(eta_n,params(4)-1,params(5)-1)';
            
            % Nodes
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount); clf;
            else
                figure('Visible', 'Off');
            end
            
            subplot(211)
            imagesc(xc_plot, zc_plot, log10(eta));
            shading flat, axis xy image, colorbar;
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            hold off
            %         title(['Viscosity (log10) at' TimeLabel])
            title(['min = ', num2str(min(eta(eta>0)), '%2.4e'), ' max = ', num2str(max(eta(eta>0)), '%2.4e')])
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1, xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minEta', 'var') caxis([minEta maxEta]); end
            
            subplot(212)
            imagesc(xc_plot, zc_plot, log10(eta_n));
            shading flat, axis xy image, colorbar;
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            hold off
            %         title(['Viscosity (log10) at' TimeLabel])
            title(['min = ', num2str(min(eta_n(eta_n>0)), '%2.4e'), ' max = ', num2str(max(eta_n(eta_n>0)), '%2.4e')])
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minEta', 'var') caxis([minEta maxEta]); end
            
            drawnow
            
            if printfig == 1
                print([path, './Fig_Viscosity/Fig_Viscosity', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % density plot
        %--------------------------------------------------
        if (rho_plot==1)
            rho_s = hdf5read(filename,'/Vertices/rho_s');  rho_s = cast(rho_s, 'double');
            rho_n = hdf5read(filename,'/Centers/rho_n');  rho_n = cast(rho_n, 'double');
            rho_n = hdf5read(filename,'/Centers/rho_n');
            rho_n = cast(rho_n, 'double');
            
            min(rho_n(:))
            max(rho_n(:))
            
            min(rho_s(:))
            max(rho_s(:))
            
            
            % Nodes
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount); clf;
            else
                figure('Visible', 'Off');
            end
            
            subplot(211)
            imagesc(xc_plot, zc_plot, (reshape(rho_s,nx,nz)'));
            shading flat, axis xy image, colorbar;
            title(['Density  nodes at' TimeLabel])
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minRho', 'var') caxis([minRho maxRho]); end
            
            subplot(212)
            imagesc(xc_plot, zc_plot, (reshape(rho_n,nx-1,nz-1)'));
            shading flat, axis xy image, colorbar;
            title(['Density  centers at' TimeLabel])
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minRho', 'var') caxis([minRho maxRho]); end
            
            if printfig == 1
                print([path, './Fig_Density/Fig_Density', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        
        
        %--------------------------------------------------
        % plot velocities errors
        %--------------------------------------------------
        if ( vel_plot == 1 )
            
            VxNodes.Vx = hdf5read(filename,'/VxNodes/Vx');
            VzNodes.Vz = hdf5read(filename,'/VzNodes/Vz');
            VxNodes.Vx = cast(VxNodes.Vx, 'double');
            VzNodes.Vz = cast(VzNodes.Vz, 'double');
            VxA = reshape(VxNodes.Vx,nx,nz+1)';
            VzA = reshape(VzNodes.Vz,params(4)+1,params(5))';
            
            
            %         figure(23)
            %         plot( xvz_plot(2:end-1), VzA(1,2:end-1));
            % %         plot( VxA(2:end-1,end), zvx_plot(2:end-1));
            %
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            clf
            % Vx
            if ( H>L )
                subplot(1,2,1)
            else
                subplot(2,1,1)
            end
            %         colormap('gray')
            
            Vx = VxA(2:end-1,:);
%             imagesc(xg_plot, zvx_plot(2:end-1), Vx*1000*3600*24*365);
            imagesc(xg_plot, zvx_plot(2:end-1), Vx);
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['Vx at' TimeLabel, ' min = ', num2str(min(Vx(:))), ' max = ', num2str(max(Vx(:)))])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minVx', 'var') caxis([minVx maxVx]); end
            
            if ( H>L )
                subplot(1,2,2)
            else
                subplot(2,1,2)
            end
            Vz = VzA(:,2:end-1);
%             imagesc(xvz_plot(2:end-1),zg_plot,Vz*1000*3600*24*365);
            imagesc(xvz_plot(2:end-1),zg_plot,Vz);

            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['Vz at' TimeLabel, ' min = ', num2str(min(Vz(:))), ' max = ', num2str(max(Vz(:)))])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minVz', 'var') caxis([minVz maxVz]); end
%             caxis([-2.5 2.5])
            
            if printfig == 1
                print([path, './Fig_Velocity/Fig_Velocity', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
            %         if print2screen == 1
            %             figCount = figCount +1;
            %             figure(figCount), clf
            %         else
            %             figure('Visible', 'Off')
            %         end
            %
            %         clf
            %
            %         dVx = hdf5read(filename,'/VxNodes/dVx');
            %         dVz = hdf5read(filename,'/VzNodes/dVz')
            %         dVx = cast(dVx, 'double');
            %         dVz = cast(dVz, 'double');
            %
            %         % Vx
            %          if ( H>L )
            %             subplot(1,2,1)
            %         else
            %             subplot(2,1,1)
            %          end
            % %         colormap('gray')
            %         Vx = reshape(dVx,nx,nz+1)';
            %         Vx = Vx(2:end-1,:);
            %         imagesc(xg_plot, zvx_plot(2:end-1), Vx);
            %         shading flat,axis xy image, colorbar;
            %         xlabel(xLabel), ylabel(zLabel);
            %         title(['dVx at' TimeLabel, ' min = ', num2str(min(Vx(:))), ' max = ', num2str(max(Vx(:)))])
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minVx', 'var') caxis([minVx maxVx]); end
            %
            %          if ( H>L )
            %             subplot(1,2,2)
            %         else
            %             subplot(2,1,2)
            %         end
            %         Vz = reshape(dVz,params(4)+1,params(5))';
            %         Vz = Vz(:,2:end-1);
            %         imagesc(xvz_plot(2:end-1),zg_plot,Vz);
            %         shading flat,axis xy image, colorbar;
            %         xlabel(xLabel), ylabel(zLabel);
            %         title(['dVz at' TimeLabel, ' min = ', num2str(min(Vz(:))), ' max = ', num2str(max(Vz(:)))])
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minVz', 'var') caxis([minVz maxVz]); end
            
            %        figure(90), hold on
            %        plot(time/1e6/3600/24/365.25,min(Vx(:))*100*365.25*24*3600,'.k')
            %        title(['normal crust min. Vx = ',num2str(min(Vx(:))*100*365.25*24*3600,'%2.2e'),'  [cm/y]'])
            %        ylabel('min. Vx [cm/y]')
            %        xlabel('time [My]')
            
            
%             figure(90), hold on
%             plot(time/1e3/3600/24/365.25,min(Vx(:))*100*365.25*24*3600,'.k')
%             title(['weak crust min. Vx = ',num2str(min(Vx(:))*100*365.25*24*3600,'%2.2e'),'  [cm/y]'])
%             ylabel('min. Vx [cm/y]')
%             xlabel('time [ky]')
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ConvTest == 1
                
                g.xmin = -1;
                g.xmax = 1;
                g.ymin = -1;
                g.ymax = 1;
                
                vxa = zeros(size(Vx));
                for i=1:nx
                    for j=1:nz-1
                        sol = eval_anal_Dani( xg_coord(i), zvx_coord(j+1) );
                        vxa(j,i) = sol.vx;
                    end
                end
                
                vza = zeros(size(Vz));
                for i=1:nx-1
                    for j=1:nz
                        sol = eval_anal_Dani( xvz_coord(i+1), zg_coord(j) );
                        vza(j,i) = sol.vz;
                    end
                end
                
                figure(10)
                subplot(1,2,1)
                imagesc(xg_coord, zvx_coord(2:end-1), (vxa))
                shading flat,axis xy image, colorbar;
                
                subplot(1,2,2)
                imagesc(xvz_coord(2:end-1), zg_coord, (vza))
                shading flat,axis xy image, colorbar;
                
                if printfig == 1
                    print([path, './Fig_Velocity/Fig_Velocity_anal', num2str(istep,'%05d'),file_suffix], '-dpng', res)
                    close all
                end
                
                ErrVx = abs(vxa-Vx);
                ErrVz = abs(vza-Vz);
                
                figCount = figCount +1;
                figure(figCount)
                subplot(1,2,1)
                imagesc(xg_coord, zvx_coord(2:end-1), (ErrVx))
                shading flat,axis xy image, colorbar;
                
                subplot(1,2,2)
                imagesc(xvz_coord(2:end-1), zg_coord, log10(abs(ErrVx)))
                shading flat,axis xy image, colorbar;
                
                nrmVx = sum(sum(abs(ErrVx(2:end-1,:)))*dx*dz)
                nrmVz = sum(sum(abs(ErrVz(:,2:end-1)))*dx*dz)
                
                if printfig == 1
                    print([path, './Fig_Velocity/Fig_Velocity_Error', num2str(istep,'%05d'),file_suffix], format, res)
                    close all
                end
                
            end
            
        end
        
        %--------------------------------------------------
        % plot velocity vectors
        %--------------------------------------------------
        if ( vel_vectors == 1 )
            
            VxNodes.Vx = hdf5read(filename,'/VxNodes/Vx'); VxNodes.Vx = cast(VxNodes.Vx, 'double');
            VzNodes.Vz = hdf5read(filename,'/VzNodes/Vz'); VzNodes.Vz = cast(VzNodes.Vz, 'double');
            Vx = reshape(VxNodes.Vx,nx,nz+1)';
            Vz = reshape(VzNodes.Vz,nx+1,nz)';
            
            % Vorticity
            Om = (Vx(2:end,:)-Vx(1:end-1,:))/dz - (Vz(:,2:end)-Vz(:,1:end-1))/dx;
            
            % Cell centered velocity
            Vx   = Vx(2:end-1,:);
            Vz   = Vz(:,2:end-1);
            Vxc  = 0.5 * (Vx(:,2:end) + Vx(:,1:end-1));
            Vzc  = 0.5 * (Vz(2:end,:) + Vz(1:end-1,:));
            
            % Stream function
            Psi1 =  -cumsum(Vzc, 2) + repmat(cumsum(Vxc(:,1)), 1, size(Vzc,2));
            Psi2 =   cumsum(Vxc, 1) - repmat(cumsum(Vzc(1,:)), size(Vzc,1), 1);
            Psi = 0.5*(Psi1 + Psi2);
            %         Psi = Psi1;
            
            %         Psi = Psi2;
            %         Psi=Psi-mean(Psi(:));    % Vx2 =  diff(Psi,1,1); % Vz2 = -diff(Psi,1,2);
            
            %         % Taneda
            %         Psi = abs(Psi); %Psi=Psi./max(Psi(:));
            Psi = Psi/100;
            %         Psi = Psi -mean(Psi(:))
            PsiContours =  [-1.4:0.04:1.2];
            PsiContours2 = [-0.70:0.015:-0.55];
            VContours    = [1.09:0.01:1.2];
            OmContours = [-4.9:0.2:4.9];
            
            %         % For the Lid Driven test of Ghia et al, 1982
            %         Psi  = Psi/100;
            %         PsiContours = [-1e-10 -1e-7 -1e-5 -1e-4 -1e-2 -3e-2 -5e-2 -7e-2 -9e-2 -0.1 -0.11 -0.115 -0.1175 ...
            %             1e-8 1e-7 1e-6 1e-5 5e-5 1e-4 2.5e-4 5e-4 1e-3 1.5e-3 3e-3];
            %         PsiContours2 = [];
            %         OmContours = [-4:6];
            
            
            %         Om = abs(Om);
            
            %         [Om]=Compute_Vorticity(nx,nz,Vx,Vz,dx,dz);
            %         [Psi]=Compute_stream(nx,nz,Om,dx,dz);
            
            dpsidx2 = diff(Psi,1,2);
            dpsidx2 = diff(dpsidx2,1,2)/dx^2;
            
            dpsidz2 = diff(Psi,1,1);
            dpsidz2 = diff(dpsidz2,1,1)/dz^2;
            
            
            Om_diff = dpsidx2(2:end-1,:) + dpsidz2(:,2:end-1);
            
            %         if crop == 1
            %             [Psi,  xc_plot, zc_plot] = CropCellArray( Psi,  xc_plot, zc_plot, lim );
            %             [Vxc,  xc_plot, zc_plot] = CropCellArray( Vxc,  xc_plot, zc_plot, lim );
            %             [Vzc,  xc_plot, zc_plot] = CropCellArray( Vzc,  xc_plot, zc_plot, lim );
            %         end
            
%             VizGrid = PhaseMap( filename, VizGrid );
%             if crop == 1
%                 [ VizGrid.ph, VizGrid.x_plot, VizGrid.z_plot ] = CropCellArray( VizGrid.ph, VizGrid.x_plot, VizGrid.z_plot, lim );
%             end

            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            eta = 0.25*(eta(1:end-1,1:end-1) + eta(2:end,1:end-1) + eta(1:end-1,2:end) + eta(2:end,2:end));
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            ezz = -exx;
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            %         sII = 2*eII.*eta;
            
            sxxd = (reshape(sxxd,params(4)-1,params(5)-1)');
            szzd = -sxxd;
            sxz  = (reshape(sxz, params(4)  ,params(5)  )');
            sII  = sqrt( 0.5*(2*sxxd.^2 + 0.5*(sxz(1:end-1,1:end-1).^2 + sxz(2:end,1:end-1).^2 + sxz(1:end-1,2:end).^2 + sxz(2:end,2:end).^2 ) ) );
            
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            HSc = 2*sxzc.*exzc + sxxd.*exx + szzd.*ezz;
            minHS=min(HSc(:));
            maxHS=max(HSc(:));
%             
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            %         if ( H>L ) subplot(121); else subplot(211); end
            Vmag = sqrt(Vxc.^2 + Vzc.^2);
            
            p=imagesc( xc_plot, zc_plot, Vmag*(365.25*3600*24*100) );
%             p=imagesc( xc_plot, zc_plot, Vxc*(365.25*3600*24*100) );
%             p=imagesc( xc_plot, zc_plot, log10(eII) );

            alpha(p,0.7);
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, Vmag*(365.25*3600*24*100), 'Visible', 'off' );
%             imagesc( xc_plot, zc_plot, log10(eII), 'Visible', 'off' );
%             imagesc( xc_plot, zc_plot, Vxc, 'Visible', 'off' );

            shading flat,axis xy image, colorbar('SouthOutside');
            
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            hv=quiver(xc_plot(1:step:end), zc_plot(1:step:end), Vxc(1:step:end,1:step:end), Vzc(1:step:end,1:step:end));
            set(hv, 'Color', 'k', 'LineWidth', 0.5);
            set(hv, 'MaxHeadSize',2.0, 'AutoScaleFactor', 2);
%             caxis([-.5 .75])
                        hold off


            xlabel(xLabel, 'FontSize', 15), ylabel(zLabel, 'FontSize', 15)
            title(['Velocity field at' TimeLabel], 'FontSize', 15)
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            set(gca, 'FontSize', 15);
%             p=imagesc( xc_plot, zc_plot, Vmag )
%             colorbar
%             alpha(p,0.7);
%             [c,hb]=contour( xc_plot, zc_plot, Vmag, VContours );
%             
%             
%             %
%             %          p=imagesc( xc_plot(2:end-1), zc_plot(2:end-1), Om_diff);
%             %          alpha(p,0.7);
%             %         [c,hb]=contour( xc_plot(2:end-1), zc_plot(2:end-1), Om_diff, OmContours );
%             
%             %         p=imagesc( xg_plot, zg_plot, Om );
%             %         alpha(p,0.7);
%             %         [c,hb]=contour( xg_plot, zg_plot, Om, OmContours );
%             %         caxis([-5 7])
%             %                 caxis([-0.05 0.05])
%             
%             %                 text_handle = clabel(c,hb);
%             % set(text_handle,'BackgroundColor',[1 1 .6],...
%             %     'Edgecolor',[.7 .7 .7])
%             %                 set(hb,'ShowText','on')
%             
%             set(hb, 'Color', 'k', 'LineWidth', 0.8);
%             shading flat,axis xy image, colorbar;
%             
%             hv=quiver(xc_plot(1:step:end), zc_plot(1:step:end), Vxc(1:step:end,1:step:end), Vzc(1:step:end,1:step:end));
%             set(hv, 'Color', 'k', 'LineWidth', 1);
%             
%             hold off
%             if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%             xlabel(xLabel), ylabel(zLabel);
%             %         title(['Velocity Magnitude at' TimeLabel])
%             title(['Vorticity at' TimeLabel])
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %         if ( H>L ) subplot(122); else subplot(212); end
            %         shading flat
            %         hold on
            %
            %         p=imagesc( xg_plot, zg_plot, Psi );
            % %         alpha(p,0.7);
            %         colorbar
            % %                 caxis([-1 1 ])
            %
            %         [c,hb]=contour( xc_plot, zc_plot, Psi, PsiContours );  %
            %         set(hb, 'Color', 'k', 'LineWidth', 0.8);
            %
            %         [c2,hb2]=contour( xc_plot, zc_plot, Psi, PsiContours2 );
            %         set(hb2, 'Color', 'w', 'LineWidth', 0.8, 'LineStyle', '-');
            %
            % %
            % %         caxis([-0.08 3e-3])
            % %
            % %         VizGrid.ph(VizGrid.ph==2) = 1;
            % %         for i = 1:1
            % %             if ( sum (sum (VizGrid.ph == i-1)) ~= 0)
            % %                 Vizbuf = VizGrid.ph;
            % %                 Vizbuf(VizGrid.ph == i) = 1;
            % %                 Vizbuf(VizGrid.ph ~= i) = 69;
            % %                 [c, hb] = contour(VizGrid.x_plot, VizGrid.z_plot, Vizbuf, [ 1 1]);
            % %                 set(hb, 'Color', 'w', 'LineWidth', 2);
            % %             end
            % %         end
            %
            %         hold off
            %         colorbar
            %         xlabel(xLabel), ylabel(zLabel);
            %         title(['Stream function at' TimeLabel])
            %         axis image;
            
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            
            if printfig == 1
                print([path, './Fig_VelocityMagnitude/Fig_VelocityMagnitude', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
        end
        
        %--------------------------------------------------
        % plot velocity divergence
        %--------------------------------------------------
        if ( vel_divergence == 1 )
            
            divV     = hdf5read(filename,'/Centers/divu');     divV    = cast(divV,    'double'); divV    = reshape(divV   , params(4)-1,params(5)-1)';
            divV_el  = hdf5read(filename,'/Centers/divu_el');  divV_el = cast(divV_el, 'double'); divV_el = reshape(divV_el, params(4)-1,params(5)-1)';
            divV_pl  = hdf5read(filename,'/Centers/divu_pl');  divV_pl = cast(divV_pl, 'double'); divV_pl = reshape(divV_pl, params(4)-1,params(5)-1)';
            divV_r   = hdf5read(filename,'/Centers/divu_r');   divV_r  = cast(divV_r,  'double'); divV_r  = reshape(divV_r , params(4)-1,params(5)-1)';
            divV_th  = hdf5read(filename,'/Centers/divu_th');  divV_th = cast(divV_th, 'double'); divV_th = reshape(divV_th, params(4)-1,params(5)-1)';
            divV_net = divV - divV_el - divV_pl - divV_r - divV_th;
             
            if crop == 1
                [divV, xc_plot, zc_plot] = CropCellArray( divV, xc_plot, zc_plot, lim );
            end
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            
            subplot(3,2,1)
            imagesc( xc_plot, zc_plot, divV_el )
            shading flat,axis xy image, colorbar;
            xlabel('X'), ylabel('Z')
            title(['Elastic divergence at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('mindiv', 'var'), caxis([mindiv maxdiv]); end
            
            subplot(3,2,2)
            imagesc( xc_plot, zc_plot, divV_pl )
            shading flat,axis xy image, colorbar;
            xlabel('X'), ylabel('Z')
            title(['Plastic divergence at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('mindiv', 'var'), caxis([mindiv maxdiv]); end
            
                        
            subplot(3,2,3)
            imagesc( xc_plot, zc_plot, divV_th )
            shading flat,axis xy image, colorbar;
            xlabel('X'), ylabel('Z')
            title(['Thermal divergence at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('mindiv', 'var'), caxis([mindiv maxdiv]); end
            
            subplot(3,2,4)
            imagesc( xc_plot, zc_plot, divV_r )
            shading flat,axis xy image, colorbar;
            xlabel('X'), ylabel('Z')
            title(['Reaction divergence at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('mindiv', 'var'), caxis([mindiv maxdiv]); end
            
            subplot(3,2,5)
            imagesc( xc_plot, zc_plot, divV )
            shading flat,axis xy image, colorbar;
            xlabel('X'), ylabel('Z')
            title(['Velocity divergence at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('mindiv', 'var'), caxis([mindiv maxdiv]); end
            
            subplot(3,2,6)
            imagesc( xc_plot, zc_plot, divV_net )
            shading flat,axis xy image, colorbar;
            xlabel('X'), ylabel('Z')
            title(['Net divergence at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('mindiv', 'var'), caxis([mindiv maxdiv]); end
            
            if printfig == 1
                print([path, './Fig_VelocityDivergence/Fig_VelocityDivergence', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
        end
        
        
        %--------------------------------------------------
        % plot pressure
        %--------------------------------------------------
        if (pre_plot==1)
            P  = hdf5read(filename,'/Centers/P');
            P  = cast(P , 'double');
            P = reshape(P,params(4)-1,params(5)-1)';
            
            %         dP  = hdf5read(filename,'/Centers/dP');
            %         dP  = cast(dP , 'double');
            %         dP = reshape(dP,params(4)-1,params(5)-1)';
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            %         P = P-mean(P(:));
            hold on
            imagesc( xc_plot, zc_plot, P)
            %         imagesc( 0.5*(xc_plot(1:end-1)+xc_plot(2:end)), zc_plot, diff(P,1,2)/(xc_plot(2)-xc_plot(1)))
            %         caxis([-10000 10000])
            shading flat,axis xy image, colorbar;
            
            %         AddCompoContours( filename, VizGrid, crop, lim  )
            %         imagesc( xc_plot, zc_plot, P, 'visible', 'off')
            %         set(hb, 'Color', 'w', 'LineWidth', 2);
            
            hold off
            %         CMAP = colormap('jet');
            %         CMAP(1,:) = [1 1 1];
            %         colormap(CMAP);
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minP', 'var') caxis([minP maxP]); end
            xlabel(xLabel), ylabel(zLabel);
            %         title(['Pressure at \gamma =' num2str(2* params(1) )])
            %                 title(['Pressure at \gamma =' num2str(2* params(1) * 365*24*3600)])
            
            %         if print2screen == 1
            %              figCount = figCount +1;
            %             figure(figCount), clf
            %         else
            %             figure('Visible', 'Off')
            %         end
            % %         P = P-mean(P(:));
            %         hold on
            %         imagesc( xc_plot, zc_plot, dP)
            %         title('dp')
            %         shading flat,axis xy image, colorbar;
            %
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         xlabel(xLabel), ylabel(zLabel);
            
            if printfig == 1
                print([path,'./Fig_Pressure/Fig_Pressure',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ConvTest == 1
                
                g.xmin = -1;
                g.xmax = 1;
                g.ymin = -1;
                g.ymax = 1;
                
                pa = zeros(size(P));
                for i=1:nx-1
                    for j=1:nz-1
                        sol = eval_anal_Dani( xc_coord(i), zc_coord(j) );
                        pa(j,i) = sol.P;
                    end
                end
                
                ErrP = P-pa;
                
                figCount = figCount +1;
                figure(figCount)
                imagesc(xc_coord,zc_coord,log10(abs(ErrP)));
                shading flat,axis xy image, colorbar;
                
                if printfig == 1
                    print([path,'./Fig_Pressure/Fig_Pressure_Error',num2str(istep,'%05d'),file_suffix], format, res)
                    close all
                end
            end
        end
        
        %--------------------------------------------------
        % plot stress and strain rate invariants
        %--------------------------------------------------
        if ( stress_inv == 1 )
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            szzd  = hdf5read(filename,'/Centers/szzd'); szzd = cast(szzd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            ezzd  = hdf5read(filename,'/Centers/ezzd'); ezzd = cast(ezzd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            eta = 0.25*(eta(1:end-1,1:end-1) + eta(2:end,1:end-1) + eta(1:end-1,2:end) + eta(2:end,2:end));
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            ezz = (reshape(ezzd,params(4)-1,params(5)-1)');
            
            exz  = (reshape(exz, params(4)  ,params(5)  )');
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            eII  = sqrt(1/2*(exx.^2+ezz.^2+(-exx-ezz).^2) + exzc.^2);
%             eII = sqrt( 0.5*( exx.^2 + ezz.^2 + (-exx-ezz).^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            %         sII = 2*eII.*eta;
            
            sxxd = (reshape(sxxd,params(4)-1,params(5)-1)');
            szzd = (reshape(szzd,params(4)-1,params(5)-1)');
            sxz  = (reshape(sxz, params(4)  ,params(5)  )');
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            sII  = sqrt(1/2*(sxxd.^2+szzd.^2+(-sxxd-szzd).^2) + sxzc.^2);

%             sII  = sqrt( 0.5*(sxxd.^2 + szzd.^2 + (-sxxd-szzd).^2 + 0.5*(sxz(1:end-1,1:end-1).^2 + sxz(2:end,1:end-1).^2 + sxz(1:end-1,2:end).^2 + sxz(2:end,2:end).^2 ) ) );
            
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            HSc = 2*sxzc.*exzc + sxxd.*exx + szzd.*ezz;
            minHS=min(HSc(:));
            maxHS=max(HSc(:));
            
            load('roma.mat')
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
                        
            % eII
            if ( H>L )
                subplot(1,2,1)
            else
                subplot(2,1,1)
            end
<<<<<<< HEAD
%             colormap(turbo);
=======
            colormap(flipud(roma));
>>>>>>> main
            imagesc( xc_plot, zc_plot, log10(eII) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, log10(eII), 'visible', 'off' )
            hold off
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['eII at' TimeLabel, ' min = ', num2str((min(eII(:))), '%2.2e'  ), ' s^{-1}', ' max = ', num2str((max(eII(:))), '%2.2e'  ), ' s^{-1}' ])
            if exist('minEii', 'var') caxis([minEii maxEii]); end
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
            
            % sII
            if ( H>L )
                subplot(1,2,2)
            else
                subplot(2,1,2)
            end
            imagesc( xc_plot, zc_plot, (sII) )
            shading flat,axis xy image, colorbar;
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (sII), 'visible', 'off' )
            hold off
            xlabel(xLabel), ylabel(zLabel);
            title(['sII at' TimeLabel, ' min = ', num2str((min(sII(:))), '%2.2e' ), ' Pa', ' max = ', num2str((max(sII(:))), '%2.2e' ), ' Pa' ])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minSii', 'var') caxis([(minSii) (maxSii)]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
            
%                         colormap('parula')

            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/Fig_Stress_StrainRate',num2str(istep,'%05d'),file_suffix], format, res)
                close all
                
            end
        end
        
        %--------------------------------------------------
        % plot stress and strain rate invariants evolution
        %--------------------------------------------------
        if ( stress_evol == 1 )
                       
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            szzd  = hdf5read(filename,'/Centers/szzd'); szzd = cast(szzd, 'double');
            P     = hdf5read(filename,'/Centers/P');    P    = cast(P, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            ezzd  = hdf5read(filename,'/Centers/ezzd'); ezzd = cast(ezzd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            eta = 0.25*(eta(1:end-1,1:end-1) + eta(2:end,1:end-1) + eta(1:end-1,2:end) + eta(2:end,2:end));
            
            exxd = reshape(exxd,params(4)-1,params(5)-1)';
            ezzd = reshape(ezzd,params(4)-1,params(5)-1)';
            exz  = reshape(exz, params(4)  ,params(5)  )';
            
            sxxd = reshape(sxxd,params(4)-1,params(5)-1)';
            P    = reshape(P,params(4)-1,params(5)-1)';
            szzd = reshape(szzd,params(4)-1,params(5)-1)';
            sxz  = reshape(sxz, params(4)  ,params(5)  )';
            
            syyd = -(sxxd + szzd);
            eyyd = -(exxd + ezzd);
            
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            HSc = 2*sxzc.*exzc + sxxd.*exxd + szzd.*ezzd;
            
            sII  = sqrt(1/2*(sxxd.^2 + szzd.^2 + syyd.^2) + sxzc.^2);
            eII  = sqrt(1/2*(exxd.^2 + ezzd.^2 + eyyd.^2) + exzc.^2);
            
            sII = sII(sII>1e-8);
            vol = length(sII)*dx*dz;
            stress = sum(sII(:).*(dx*dz))/vol;
            press  = sum(P(:).*(dx*dz))/vol;

            if istep == 0
                stress = 0;
            end
            figure(90), 
            if istep==0 || icount==1, clf; end
            subplot(211), hold on
            plot(time/1e3/3600/365/24, (stress), 'k.'), ylabel('Tii')
            subplot(212), hold on
            plot(time/1e3/3600/365/24, (press), 'k.'), ylabel('P')
            
            short(icount)  =  strain;
            siivec(icount) = sum(sII(:).*(dx*dz))/vol;
            eiivec(icount) = sum(eII(:).*(dx*dz))/vol;
            eiimaxvec(icount) = max(eII(:));
            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/Fig_Stress_StrainRate',num2str(istep,'%05d'),file_suffix], format, res)
                close all
                
            end
            
            if istep == iend
                save('MaxStressStrainRate', 'siivec', 'eiivec', 'eiimaxvec', 'short')
            end
            
        end
        
        %--------------------------------------------------
        % plot stress and strain rate invariants evolution
        %--------------------------------------------------
        if ( temp_evol == 1 )
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            
            T  = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = (reshape(T,params(4)-1,params(5)-1)');
            
            VizGrid = PhaseMap( filename, VizGrid );
            
            Tpost = T;
            Tpost(VizGrid.ph~=7 & VizGrid.ph~=8) = 0;
            Tpre = T;
            Tpre(VizGrid.ph~=1 & VizGrid.ph~=2) = 0;
            T_postrift(icount) = max(Tpost(:));
            T_prerift(icount) = max(Tpre(:));
            

            if print2screen == 1
            figure(90)
            hold on            
            plot(time/1e6/3600/365/24, T_prerift(icount)-273.15, '.r')
            plot(time/1e6/3600/365/24, T_postrift(icount)-273.15, '.b')
            
            figure(91), clf
            imagesc( xc_plot, zc_plot,Tpre-273.15 )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, T_prerift-273.15, 'visible', 'off' )
            hold off
            shading flat,axis xy image, colorbar;
            caxis([0 900])
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            end
            short(icount)  =  strain;
            timevec(icount) = time/1e6/3600/365/24;
         
            if istep == iend && printfig==1
                save('MaxTemperatureSediments', 'timevec', 'T_prerift', 'T_postrift', 'short')
            end
            
        end
        
        %--------------------------------------------------
        % plot stress components
        %--------------------------------------------------
        if ( stress_plot == 1 )
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            szzd  = hdf5read(filename,'/Centers/szzd'); szzd = cast(szzd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            
            sxxd = (reshape(sxxd,params(4)-1,params(5)-1)');
            szzd = (reshape(szzd,params(4)-1,params(5)-1)');
            sxz  = (reshape(sxz, params(4)  ,params(5)  )');
            P  = hdf5read(filename,'/Centers/P');
            P  = cast(P , 'double');
            P = reshape(P,params(4)-1,params(5)-1)';
            mean(P(:))
            %         P = P-mean(P(:));
            %         Sxx = sxxd - P;
            %         Szz = szzd - P;
            
            %         figure(90), hold on
            %         plot(time(icount), s_num_xx(icount), 'dg')
            %         plot(time(icount), s_num_zz(icount), 'sb')
            %         plot(time(icount), s_anal(icount), 'or')
            
            if ConvTest == 1
                
                g.xmin = -1;
                g.xmax =  1;
                g.ymin = -1;
                g.ymax =  1;
                
                Sxxa = zeros(size(P));
                Szza = zeros(size(P));
                for i=1:nx-1
                    for j=1:nz-1
                        sol = eval_anal_Dani( xc_coord(i), zc_coord(j) );
                        Sxxa(j,i) = sol.sxx;
                        Szza(j,i) = sol.szz;
                    end
                end
                
                ErrSxx = Sxx-Sxxa;
                ErrSzz = Szz-Szza;
                
                figure(98)
                imagesc(xc_coord,zc_coord,log10(abs(Sxxa-Sxx)));
                %                         imagesc(xc_coord,zc_coord,Sxxa);
                
                shading flat,axis xy image, colorbar;
                
                if printfig == 1
                    print([path,'./Fig_StressComponents/Fig_Sxx_Error',num2str(istep,'%05d'),file_suffix], format, res)
                    close all
                end
                nrmSxx = sum(sum(abs(ErrSxx))*dx*dz)
                %             nrmSzz = sum(sum(abs(ErrSzz))*dx*dz)
            end
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            % exxd
            if ( H>L )
                subplot(2,3,1)
            else
                subplot(3,2,1)
            end
            hold on
            %         sxxd(sxxd<0)=0;
            imagesc(xc_plot, zc_plot, sxxd)
%             contour(xc_plot, zc_plot, sxxd, [.0 .0], 'k')
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['\tau_{xx} at' TimeLabel])
            hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minSxx', 'var') caxis([minTxx maxTxx]); end
            
            % ezzd
            if ( H>L )
                subplot(2,3,2)
            else
                subplot(3,2,2)
            end
            hold on
            imagesc(xc_plot, zc_plot, szzd)
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['\tau{zz} at' TimeLabel])
            hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minSxx', 'var') caxis([-maxTxx -minTxx]); end
            
            % exz
            if ( H>L )
                subplot(2,3,3)
            else
                subplot(3,2,3)
            end
            hold on
            %         imagesc(xg_plot, zg_plot, log10(abs(sxz)));
            imagesc(xg_plot, zg_plot, sxz);
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['\tau{xz} at' TimeLabel])
            hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minSxz', 'var') caxis([minSxz maxSxz]); end
            
            % exxd
            if ( H>L )
                subplot(2,3,4)
            else
                subplot(3,2,4)
            end
            hold on
            imagesc(xc_plot, zc_plot, sxxd-P)
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['\sigma_{xx} at' TimeLabel])
            hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minSxx', 'var') caxis([minSxx maxSxx]); end
            
            % ezzd
            if ( H>L )
                subplot(2,3,5)
            else
                subplot(3,2,5)
            end
            hold on
            imagesc(xc_plot, zc_plot, -sxxd-P)
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['\sigma_{zz} at' TimeLabel])
            hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minSxx', 'var') caxis([-maxSxx -minSxx]); end
            
            % exz
            if ( H>L )
                subplot(2,3,6)
            else
                subplot(3,2,6)
            end
            hold on
            imagesc(xc_plot, zc_plot, P);
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['P at' TimeLabel])
            hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minSxz', 'var') caxis([minSxz maxSxz]); end
            
            
            if printfig == 1
                print([path,'./Fig_StressComponents/Fig_Stress_components',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot strain rate components
        %--------------------------------------------------
        if ( srate_plot == 1 )
            
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            ezzd  = hdf5read(filename,'/Centers/ezzd'); ezzd = cast(ezzd, 'double');
            %         ezzd  = hdf5read(filename,'/Centers/ezzd'); ezzd = cast(ezzd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            
            exxd = (reshape(exxd,params(4)-1,params(5)-1)');
            ezzd = (reshape(ezzd,params(4)-1,params(5)-1)');
            exz  = (reshape(exz, params(4)  ,params(5)  )');
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            % exxd
            if ( H>L )
                subplot(1,3,1)
            else
                subplot(3,1,1)
            end
            hold on
            imagesc(xc_plot, zc_plot, log10(abs(exxd)))
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['\epsilon_{xx} at' TimeLabel])
            hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minExx', 'var') caxis([minExx maxExx]); end
            
            % ezzd
            if ( H>L )
                subplot(1,3,2)
            else
                subplot(3,1,2)
            end
            hold on
            imagesc(xc_plot, zc_plot, log10(abs(ezzd)))
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['\epsilon_{zz} at' TimeLabel])
            hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minExx', 'var') caxis([minExx maxExx]); end
            
            % exz
            if ( H>L )
                subplot(1,3,3)
            else
                subplot(3,1,3)
            end
            hold on
            imagesc(xg_plot, zg_plot, log10(abs(exz)));
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['\epsilon_{xz} at' TimeLabel])
            hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minExz', 'var') caxis([minExz maxExz]); end
            
            if printfig == 1
                print([path,'./Fig_StrainRateComponents/Fig_Srate_components',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot phases on visualisation grid
        %--------------------------------------------------
        if(phase_on_grid==1)
            
            VizGrid = PhaseMap( filename, VizGrid );
            
            %         % Ball fraction
            %         frac_boules   = sum(sum(VizGrid.ph==2)) / (size(VizGrid.ph,1)*size(VizGrid.ph,2));
            %         mass_err      = abs(frac_boules - 0.3) / 0.3 * 100;
            %         error(icount) = mass_err;
            %         gamma(icount) = 2* params(1) / 1;
            
            if crop == 1
                [ VizGrid.ph, VizGrid.x_plot, VizGrid.z_plot ] = CropCellArray( VizGrid.ph, VizGrid.x_plot, VizGrid.z_plot, lim );
            end
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            if strcmp(color_map,'Benoit_MLPS')
                jet = colormap(map);
                ind = 63;
                jet(ind:end,:) = repmat([1 1 1], size(jet,1)-ind+1,1);
                colormap(jet)
                p = imagesc(VizGrid.x_plot, VizGrid.z_plot, abs(VizGrid.ph-max(VizGrid.ph(:))));
            else
                p = imagesc(VizGrid.x_plot, VizGrid.z_plot, VizGrid.ph);
            end
            xlabel(xLabel), ylabel(zLabel)
            title(['Phases at' TimeLabel]);
            shading flat, axis xy image
            hold off
            colorbar;
            colormap(map)
            if strcmp(color_map,'geoffroy_ast'), caxis ([-1 size(map,1)-2]); end
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minPhase', 'var') caxis([minPhase maxPhase]); end
            
            %           set(gca, 'XTick', []);
            %           set(gca, 'YTick', []);
            %         grid(gca, 'minor')
            %         alpha (p,0.7)
            
            if printfig == 1
                if crop == 1
                    print([path, './Fig_Phases/Crop_PhasesGrid', num2str(istep,'%05d'),file_suffix], format, res)
                else
                    print([path, './Fig_Phases/Fig_PhasesGrid', num2str(istep,'%05d'),file_suffix], format, res)
                end
                close all
            end
            
            clear VizGrid.ph VizGrid.x VizGrid.z
            
        end
        
        
        %--------------------------------------------------
        % plot hydrodynamic pressure
        %--------------------------------------------------
        if (dyna_pre==1)
            
            P  = hdf5read(filename,'/Centers/P');
            P  = cast(P , 'double');
            P = reshape(P,params(4)-1,params(5)-1)';
            %         P = P - mean(P(:));
            
            rho_n = hdf5read(filename,'/Centers/rho_n');
            rho_n = cast(rho_n, 'double');
            rhoc = reshape(rho_n,params(4)-1,params(5)-1)';
            
            % Build lithostatic pressure
            dz = abs(zc_coord(2)-zc_coord(1));
            %         rhoc = 0.25*(rho(1:end-1,1:end-1) + rho(1:end-1, 2:end) + rho(2:end, 1:end-1) + rho(2:end,2:end));
            Plith = P;
            Plith(ncz, :) = rhoc(ncz,:)*gz*dz;
            for o=ncz-1:-1:1
                Plith(o, :) = Plith(o+1, :) - rhoc(o,:)*gz*dz;
            end
            Pdyn = P-Plith;
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            imagesc( xc_plot, zc_plot, Pdyn )
            shading flat,axis xy image, colorbar;
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, Pdyn, 'Visible', 'off' )
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Dynamic pressure at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minPdyn', 'var') caxis([minPdyn maxPdyn]); end
            
            if printfig == 1
                print([path,'./Fig_DynamicPressure/Fig_DynamicPressure',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end

            
%             if print2screen == 1
%                 figCount = figCount +1;
%                 figure(figCount), clf
%             else
%                 figure('Visible', 'Off')
%             end
%             
%             imagesc( xc_plot, zc_plot, Plith )
%             shading flat,axis xy image, colorbar;
%             hold on
%             if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
%             imagesc( xc_plot, zc_plot, Pdyn, 'Visible', 'off' )
%             hold off
%             xlabel(xLabel), ylabel(zLabel)
%             title(['Lithostatic pressure at' TimeLabel])
%             if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%             
            
%             if print2screen == 1
%                 figCount = figCount +1;
%                 figure(figCount), clf
%             else
%                 figure('Visible', 'Off')
%             end
% 
%             imagesc( xc_plot, 0.5*(zc_plot(2:end)+zc_plot(1:end-1)), diff(Plith,1,1)/dz*1e-5/1000*1000 )
% 
%             hold on
%             %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
%             imagesc( xc_plot, 0.5*(zc_plot(2:end)+zc_plot(1:end-1)), diff(Plith,1,1)/dz*1e-5/1000*1000, 'Visible', 'off' )
%             hold off
%             xlabel(xLabel), ylabel(zLabel)
%             title(['dPlith/dz [Pa/m] at' TimeLabel])
%             shading flat,axis xy image, colorbar; caxis([[-0.2757 -0.2747]])
%             if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%             if printfig == 1
%                 print([path,'./Fig_DynamicPressure/Fig_LithPressure',num2str(istep,'%05d'),file_suffix], format, res)
%                 close all
%             end
            
            
        end
        
        %--------------------------------------------------
        % plot additive effective strain rates
        %--------------------------------------------------
        if (srate_add==1)
           
            eII_el  = hdf5read(filename,'/Centers/eII_el');
            eII_el  = cast(eII_el , 'double');
            eII_el = reshape(eII_el,params(4)-1,params(5)-1)';
            
            eII_pl  = hdf5read(filename,'/Centers/eII_pl');
            eII_pl  = cast(eII_pl , 'double');
            eII_pl = reshape(eII_pl,params(4)-1,params(5)-1)';
            
            eII_pwl  = hdf5read(filename,'/Centers/eII_pwl');
            eII_pwl  = cast(eII_pwl , 'double');
            eII_pwl = reshape(eII_pwl,params(4)-1,params(5)-1)';
            
            eII_exp  = hdf5read(filename,'/Centers/eII_exp');
            eII_exp  = cast(eII_exp , 'double');
            eII_exp = reshape(eII_exp,params(4)-1,params(5)-1)';
            
            eII_lin  = hdf5read(filename,'/Centers/eII_lin');
            eII_lin  = cast(eII_lin , 'double');
            eII_lin = reshape(eII_lin,params(4)-1,params(5)-1)';
            
            eII_gbs  = hdf5read(filename,'/Centers/eII_gbs');
            eII_gbs  = cast(eII_gbs , 'double');
            eII_gbs = reshape(eII_gbs,params(4)-1,params(5)-1)';
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end            
            
%             subplot(311)
%             imagesc( xc_plot, zc_plot, log10(eII_el) )
%             hold on
%             if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
%             imagesc( xc_plot, zc_plot, log10(eII_el), 'Visible', 'off' )
%             shading flat,axis xy image, colorbar;
%             if exist('minEii', 'var') caxis([minEii maxEii]); end
%             shading flat,axis xy image, colorbar;
%             xlabel(xLabel), ylabel(zLabel)
%             title(['Elastic strain rate at' TimeLabel])
%             if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%             if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
%             axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
%             
%             
%             subplot(312)
%             imagesc( xc_plot, zc_plot, log10(eII_pwl) )
%             hold on
%             if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
%             imagesc( xc_plot, zc_plot, log10(eII_pwl), 'Visible', 'off' )
%             shading flat,axis xy image, colorbar;
%             if ColorFabio==1,             load('roma.mat')
%             colormap(flipud(roma));
%             end
%             if exist('minEii', 'var') caxis([minEii maxEii]); end
%             hold off
%             xlabel(xLabel), ylabel(zLabel)
%             title(['Power law strain rate at' TimeLabel])
%             if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%             if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
%             axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
%             
%             subplot(313)
%             imagesc( xc_plot, zc_plot, log10(eII_pl) )
%             hold on
%             if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
%             imagesc( xc_plot, zc_plot, log10(eII_pl), 'Visible', 'off' )
%             shading flat,axis xy image, colorbar;
%             if exist('minEii', 'var') caxis([minEii maxEii]); end
%             hold off
%             xlabel(xLabel), ylabel(zLabel)
%             title(['Plastic strain rate at' TimeLabel])
%             if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%             if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
%             axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
            
            colormap('jet');
            
            subplot(321)
            imagesc( xc_plot, zc_plot, log10(eII_el) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, log10(eII_el), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([minEii maxEii]); end
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel)
            title(['Elastic strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(323)
            imagesc( xc_plot, zc_plot, log10(eII_pl) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, log10(eII_pl), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([minEii maxEii]); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Plastic strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(322)
            imagesc( xc_plot, zc_plot, log10(eII_pwl) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, log10(eII_pwl), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([minEii maxEii]); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Power law strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
%             caxis([-34 -20])
            
            
            
            subplot(324)
            imagesc( xc_plot, zc_plot, (log10(eII_exp)) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (log10(eII_exp)), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([minEii maxEii]); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Exponential strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(325)
            imagesc( xc_plot, zc_plot, (log10(eII_lin)) )
            %         imagesc( xc_plot, zc_plot, eII_lin./eII )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (log10(eII_lin)), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([minEii maxEii]); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Linear strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
%             caxis([-26.5 -20])
            
            
            subplot(326)
            imagesc( xc_plot, zc_plot, (log10(eII_gbs)) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (log10(eII_gbs)), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([minEii maxEii]); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['GBS strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/StrainRatesAdd',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot accumulated additive strains
        %--------------------------------------------------
        if (acc_strain_add==1)
            
            eII_el  = hdf5read(filename,'/Centers/strain_el');
            eII_el  = cast(eII_el , 'double');
            eII_el = reshape(eII_el,params(4)-1,params(5)-1)';
            
            eII_pl  = hdf5read(filename,'/Centers/strain_pl');
            eII_pl  = cast(eII_pl , 'double');
            eII_pl = reshape(eII_pl,params(4)-1,params(5)-1)';
            
            eII_pwl  = hdf5read(filename,'/Centers/strain_pwl');
            eII_pwl  = cast(eII_pwl , 'double');
            eII_pwl = reshape(eII_pwl,params(4)-1,params(5)-1)';
            
            eII_lin = zeros(size(eII_pwl));
            eII_lin  = hdf5read(filename,'/Centers/strain_lin');
            eII_lin  = cast(eII_lin , 'double');
            eII_lin = reshape(eII_lin,params(4)-1,params(5)-1)';
            
            eII_exp  = hdf5read(filename,'/Centers/strain_exp');
            eII_exp  = cast(eII_exp , 'double');
            eII_exp = reshape(eII_exp,params(4)-1,params(5)-1)';
            
            eII  = hdf5read(filename,'/Centers/strain');
            eII  = cast(eII , 'double');
            eII = reshape(eII,params(4)-1,params(5)-1)';
            strmin = min(eII(:));
            strmax = max(eII(:))/4;
            
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            subplot(321)
            imagesc( xc_plot, zc_plot, log10(eII_el) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (eII_el), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([strmin strmax]); end
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel)
            title(['Elastic strain at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(323)
            imagesc( xc_plot, zc_plot, (eII_pl) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (eII_pl), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([strmin strmax]); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Plastic strain at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(322)
            imagesc( xc_plot, zc_plot, (eII_pwl) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (eII_pwl), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([strmin strmax]); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Power law strain at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(324)
            imagesc( xc_plot, zc_plot, ((eII_exp)) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, ((eII_exp)), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([strmin strmin]); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Exponential strain at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(325)
            imagesc( xc_plot, zc_plot, ((eII_lin)) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, ((eII_lin)), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([strmin strmax]); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Linear strain at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(326)
            imagesc( xc_plot, zc_plot, ((eII)) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, ((eII)), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var') caxis([strmin strmax]); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Total strain at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            %         subplot(221)
            %         imagesc( xc_plot, zc_plot, (eII_el) )
            %         shading flat,axis xy image, colorbar;
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         imagesc( xc_plot, zc_plot, (eII_el), 'Visible', 'off' )
            %         hold off
            %         xlabel(xLabel), ylabel(zLabel)
            %         title(['Acc. Elastic strain at' TimeLabel])
            %         if exist('strmin', 'var') caxis([(strmin) (strmax)]); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %
            %         subplot(223)
            %         imagesc( xc_plot, zc_plot, (abs(eII_pl)) )
            %         shading flat,axis xy image, colorbar;
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         imagesc( xc_plot, zc_plot, (abs(eII_pl)), 'Visible', 'off' )
            %         hold off
            %         xlabel(xLabel), ylabel(zLabel)
            %         title(['Acc. Plastic strain at' TimeLabel])
            %         if exist('strmin', 'var') caxis([(strmin) (strmax)]); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %
            %         subplot(222)
            %         imagesc( xc_plot, zc_plot, (abs(eII_pwl)) )
            %         shading flat,axis xy image, colorbar;
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         imagesc( xc_plot, zc_plot, (abs(eII_pwl)), 'Visible', 'off' )
            %         hold off
            %         xlabel(xLabel), ylabel(zLabel)
            %         title(['Acc. Power law strain at' TimeLabel])
            %         if exist('strmin', 'var') caxis([(strmin) (strmax)]); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %
            %         subplot(224)
            %         imagesc( xc_plot, zc_plot, ((eII_exp)) )
            %         shading flat,axis xy image, colorbar;
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         imagesc( xc_plot, zc_plot, ((eII_exp)), 'Visible', 'off' )
            %         hold off
            %         xlabel(xLabel), ylabel(zLabel)
            %         title(['Acc. Exp. strain at' TimeLabel])
            %         if exist('strmin', 'var') caxis([(strmin) (strmax)]); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            if printfig == 1
                print([path,'./Fig_AccStrain/Fig_AccStrainAdd',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot accumulated strain
        %--------------------------------------------------
        if (acc_strain==1)
            
            strain = hdf5read(filename,'/Centers/strain'); strain = cast(strain, 'double');
            strain = reshape(strain,params(4)-1,params(5)-1)';
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            
            %         jet1 = colormap('jet');
            %         ind = 1;
            %         jet1(1:ind,:) = repmat([1 1 1], ind,1);
            %         colormap(jet1)
            
            minstrain = -4;
            logstrain = log10(strain);
            %         logstrain(logstrain < 0.9*minstrain & logstrain~=-Inf ) = 0.9*minstrain;
            %         logstrain(logstrain == -Inf ) = -10;
            
            p=imagesc( xc_plot, zc_plot, logstrain );
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (log10(strain)), 'Visible', 'off' );
            shading flat,axis xy image, colorbar;
            hold off
            if exist('minStr', 'var') caxis([minStr maxStr]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end

            xlabel(xLabel, 'FontSize', 15), ylabel(zLabel, 'FontSize', 15)
            title(['Acc. strain at' TimeLabel], 'FontSize', 15)
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            set(gca, 'FontSize', 15);

            if printfig == 1
                print([path, './Fig_AccStrain/Fig_AccStrain', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot temperature
        %--------------------------------------------------
        if (temperature==1)
            
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)' - 273.15;
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            imagesc( xc_plot, zc_plot, T )
            shading flat,axis xy image, colorbar;
            colormap('jet')
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         imagesc( xc_plot, zc_plot, T, 'Visible', 'off' )
            
            imagesc( xc_plot, zc_plot, T, 'Visible', 'off' )
            hold off
            
            xlabel(xLabel), ylabel(zLabel)
            title(['T at' TimeLabel ' min = ' num2str(min(T(:))) ' C max = ' num2str(max(T(:))) ' C'])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minT', 'var') caxis([minT maxT]); end
            %         caxis([273 1000])
            
            if printfig == 1
                print([path, './Fig_Temperature/Fig_Temperature', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot for Stefan
        %--------------------------------------------------
        if (Stefan==1)
            
            show = 1;
            
            dt = time-time0; % seconds
            
            %         if istep>0
            %
            %         end
            
            % Temperature
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)';
            
            % Viscosity
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            
            % Strain rate ariant II
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*(2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ));
            
            %         sII  = hdf5read(filename,'/Centers/sII'); sII = cast(sII, 'double');
            %         eII  = hdf5read(filename,'/Centers/eII'); eII = cast(eII, 'double');
            
            %         sII = (reshape(sII,params(4)-1,params(5)-1)');
            %         eII = (reshape(eII,params(4)-1,params(5)-1)');
            
            % Total stress
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            P     = hdf5read(filename,'/Centers/P');    P    = cast(P, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            sxz = (reshape(sxz, params(4)  ,params(5)  )');
            
            sxxd = (reshape(sxxd,params(4)-1,params(5)-1)');
            P    = (reshape(   P,params(4)-1,params(5)-1)');
            sII = sqrt( 0.5*(2*sxxd.^2 + 0.5*(sxz(1:end-1,1:end-1).^2 + sxz(2:end,1:end-1).^2 + sxz(1:end-1,2:end).^2 + sxz(2:end,2:end).^2 ) ));
            
            sxx   =  sxxd - P;
            szz   = -sxxd - P;
            sxz  = (reshape(sxz, params(4)  ,params(5)  )');
            sII  = sqrt( 0.5*(2*sxxd.^2 + 0.5*(sxz(1:end-1,1:end-1).^2 + sxz(2:end,1:end-1).^2 + sxz(1:end-1,2:end).^2 + sxz(2:end,2:end).^2 ) ) );
            
            % Velocities (cell centers)
            VxNodes.Vx = hdf5read(filename,'/VxNodes/Vx'); VxNodes.Vx = cast(VxNodes.Vx, 'double');
            VzNodes.Vz = hdf5read(filename,'/VzNodes/Vz'); VzNodes.Vz = cast(VzNodes.Vz, 'double');
            Vx = reshape(VxNodes.Vx,nx,nz+1)';
            Vz = reshape(VzNodes.Vz,nx+1,nz)';
            Vx = Vx(2:end-1,:); Vx = 0.5 * (Vx(:,2:end) + Vx(:,1:end-1));
            Vz = Vz(:,2:end-1); Vz = 0.5 * (Vz(2:end,:) + Vz(1:end-1,:));
            
            AB = 1;
            
            if icount > 1
                
                file0 = [path,'Output',num2str(istep-ijump,'%05d'),file_suffix,'.gzip.h5'];
                T0  = hdf5read(file0,'/Centers/T');
                T0 = cast(T0, 'double');
                T0 = reshape(T0,params(4)-1,params(5)-1)';
                
                Wdot = sum( dz*(  abs(sxx(:,1).*Vx(:,1)) + abs(sxx(:,end).*Vx(:,end)) )) + sum(dx*(  abs(szz(1,:).*Vz(1,:)) +  abs(szz(end,:).*Vz(end,:))  ));
                fW_0 = dt*Wdot;
                dT   = T - T0;
                fE_0 = sum(sum(dx*dz*1050*2700.*abs(dT)));
                %             figure(11)%, hold on
                %             plot( strain, (Es), '.k', strain, (Ws), '+k')
                
                if icount > 2
                    
                    file_1 = [path,'Output',num2str(istep-2*ijump,'%05d'),file_suffix,'.gzip.h5'];
                    T_1    = hdf5read(file_1,'/Centers/T'); T_1    = cast(T_1, 'double');
                    T_1    = reshape(T_1,params(4)-1,params(5)-1)';
                    xg_1   = hdf5read(file_1,'/Model/xg_coord');  xg_coord  = cast(xg_coord, 'double'); dx_1 = xg_1(2) - xg_1(1);
                    zg_1   = hdf5read(file_1,'/Model/zg_coord');  zg_coord  = cast(zg_coord, 'double'); dz_1 = zg_1(2) - zg_1(1);
                    
                    dT0     = T0 - T_1;
                    fE_1    = sum(sum(dx_1*dz_1*1050*2700.*abs(dT0)));
                    
                    % Total stress
                    sxxd_1 = hdf5read(file_1,'/Centers/sxxd'); sxxd_1 = cast(sxxd_1, 'double');
                    P_1    = hdf5read(file_1,'/Centers/P');    P_1    = cast(P_1, 'double');
                    sxxd_1 = (reshape(sxxd_1,params(4)-1,params(5)-1)');
                    P_1    = (reshape(   P_1,params(4)-1,params(5)-1)');
                    sxx_1  =  sxxd_1 - P_1;
                    szz_1  = -sxxd_1 - P_1;
                    
                    % Velocities (cell centers)
                    Vx_1 = hdf5read(file_1,'/VxNodes/Vx'); Vx_1 = cast(Vx_1, 'double');
                    Vz_1 = hdf5read(file_1,'/VzNodes/Vz'); Vz_1 = cast(Vz_1, 'double');
                    Vx_1 = reshape(Vx_1,nx,nz+1)';
                    Vz_1 = reshape(Vz_1,nx+1,nz)';
                    Vx_1 = Vx_1(2:end-1,:); Vx_1 = 0.5 * (Vx_1(:,2:end) + Vx_1(:,1:end-1));
                    Vz_1 = Vz_1(:,2:end-1); Vz_1 = 0.5 * (Vz_1(2:end,:) + Vz_1(1:end-1,:));
                    Wdot_1 = sum( dz_1*(  abs(sxx_1(:,1).*Vx_1(:,1)) + abs(sxx_1(:,end).*Vx_1(:,end)) )) + sum(dx_1*(  abs(szz_1(1,:).*Vz_1(1,:)) +  abs(szz_1(end,:).*Vz_1(end,:))  ));
                    dt0    = tits.Time(icount-1) - tits.Time(icount-2);
                    fW_1   = dt0*Wdot_1;
                    
                    
                    if icount > 3
                        
                        file_2 = [path,'Output',num2str(istep-3*ijump,'%05d'),file_suffix,'.gzip.h5'];
                        T_2    = hdf5read(file_2,'/Centers/T'); T_2    = cast(T_2, 'double');
                        T_2    = reshape(T_2,params(4)-1,params(5)-1)';
                        xg_2   = hdf5read(file_2,'/Model/xg_coord');  xg_coord  = cast(xg_coord, 'double'); dx_2 = xg_2(2) - xg_2(1);
                        zg_2   = hdf5read(file_2,'/Model/zg_coord');  zg_coord  = cast(zg_coord, 'double'); dz_2 = zg_2(2) - zg_2(1);
                        
                        dT_1   = T_1 - T_2;
                        fE_2   = sum(sum(dx_2*dz_2*1050*2700.*abs(dT_1)));
                        
                        % Total stress
                        sxxd_2  = hdf5read(file_2,'/Centers/sxxd'); sxxd_2 = cast(sxxd_2, 'double');
                        P_2     = hdf5read(file_2,'/Centers/P');    P_2    = cast(P_2, 'double');
                        sxxd_2 = (reshape(sxxd_2,params(4)-1,params(5)-1)');
                        P_2    = (reshape(   P_2,params(4)-1,params(5)-1)');
                        sxx_2   =  sxxd_2 - P_2;
                        szz_2   = -sxxd_2 - P_2;
                        
                        % Velocities (cell centers)
                        Vx_2 = hdf5read(file_2,'/VxNodes/Vx'); Vx_2 = cast(Vx_2, 'double');
                        Vz_2 = hdf5read(file_2,'/VzNodes/Vz'); Vz_2 = cast(Vz_2, 'double');
                        Vx_2 = reshape(Vx_2,nx,nz+1)';
                        Vz_2 = reshape(Vz_2,nx+1,nz)';
                        Vx_2 = Vx_2(2:end-1,:); Vx_2 = 0.5 * (Vx_2(:,2:end) + Vx_2(:,1:end-1));
                        Vz_2 = Vz_2(:,2:end-1); Vz_2 = 0.5 * (Vz_2(2:end,:) + Vz_2(1:end-1,:));
                        Wdot_2    = sum( dz_2*(  abs(sxx_2(:,1).*Vx_2(:,1)) + abs(sxx_2(:,end).*Vx_2(:,end)) )) + sum(dx_2*(  abs(szz_2(1,:).*Vz_2(1,:)) +  abs(szz_2(end,:).*Vz_2(end,:))  ));
                        dt_1 = tits.Time(icount-2) - tits.Time(icount-3);
                        fW_2 = dt_1*Wdot_2;
                        
                    end
                end
            end
            
            disp('------------------')
            icount
            
            if icount == 1
                updateE = 0;
                updateW = 0;
            end
            
            if icount > 1
                updateE = fE_0;
                updateW = fW_0;
            end
            
            
            if (icount > 2 && AB>1) || (icount == 3 && AB>2)
                updateE = 3/2*fE_1 - 1/2*fE_0;
                updateW = 3/2*fW_1 - 1/2*fW_0;
            end
            %
            if icount > 3 && AB>2
                updateE = 23/12*fE_2 - 16/12*fE_1 + 5*fE_0;
                updateW = 23/12*fW_2 - 16/12*fW_1 + 5*fW_0;
            end
            
            Es   = Es + updateE;
            Ws   = Ws + updateW;
            
            if istep == 0
                
                %             figure(11), hold on
                %             stef=load('/Users/tduretz/Desktop/Cooling/Debug/ShearHeat_Data');
                %             plot(stef.Sho,(stef.E), '-k');
                %             plot(stef.Sho,(stef.W), '-.k');
                
            end
            
            % Save in arrays
            tits.W(icount)    = Ws;
            tits.E(icount)    = Es;
            tits.Time(icount) = time;
            tits.Sho(icount)  = strain;
            
            if show
                
                if print2screen == 1
                    figCount = figCount +1;
                    figure(figCount), clf
                else
                    figure('Visible', 'Off')
                end
                
                % eta
                subplot(3,1,1),
                imagesc(xc_plot,zc_plot, (sII));
                if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
                %         imagesc(xc_plot,zc_plot, (sII), 'Visible', 'off');
                %         drawnow
                if exist('minSii', 'var') caxis([minSii maxSii]); end
                shading flat;
                hold on;
                axis xy image, colorbar;
                xlabel(xLabel), ylabel(zLabel);
                title(['sII at' TimeLabel  ' min = ' num2str(min((sII(:))), '%2.2e') ' Pa max = ' num2str(max((sII(:))), '%2.2e') ' Pa'])
                if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
                %         if exist('minE', 'var') caxis([minE maxE]); end
                %         imagesc(xg_plot, zg_plot, log10(eta));
                %         shading flat
                %         title(['Viscosity  nodes (log10) at' TimeLabel  ' min = ' num2str(min((eta(:))), '%2.2e') ' Pa.s max = ' num2str(max((eta(:))), '%2.2e') ' Pa.s'])
                %         xlabel(xLabel), ylabel(zLabel);
                %         axis xy image, colorbar;
                %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
                %         if exist('minEta', 'var') caxis([minEta maxEta]); end
                
                % eII
                subplot(3,1,2),
                imagesc(xc_plot,zc_plot,log10(eII));
                %          if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
                
                %         CMAP = colormap('jet');
                %         CMAP(1,:) = [0.99 0.99 0.99];
                %         colormap(CMAP);
                %         if exist('minEii', 'var') caxis([minEii maxEii]); end
                shading flat;
                hold on;
                %         [c, hb] = contour(xc_plot,zc_plot,((eII)-eps_val>0), 1);
                %         set(hb, 'Color', 'w', 'LineWidth', 2);
                axis xy image, colorbar;
                xlabel(xLabel), ylabel(zLabel);
                title(['eII at' TimeLabel  ' min = ' num2str(min((eII(:)))) ' /s max = ' num2str(max((eII(:)))) ' /s'])
                if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
                if exist('minEii', 'var') caxis([minEii maxEii]); end
                
                subplot(3,1,3)
                imagesc(xc_plot, zc_plot, T-273);
                %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
                
                shading flat
                colorbar;%('SouthOutside')
                axis xy image;
                xlabel(xLabel), ylabel(zLabel)
                title(['T at' TimeLabel ' min = ' num2str(min(T(:))-273) ' C max = ' num2str(max(T(:))-273) ' C'])
                if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
                if exist('minT', 'var') caxis([minT maxT]); end
                
                if printfig == 1
                    print([path, './Fig_Temperature/Fig_Stefan', num2str(istep,'%05d'),file_suffix], format, res)
                    close all
                end
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %         if (istep==iend)
            %             figCount = figCount +1;
            %             figure(figCount), hold on
            %             %             stef=load('/Users/tduretz/Desktop/Cooling/Debug/ShearHeat_Data');
            %             %             plot(stef.Sho,(stef.E), '-k');
            %             %             plot(stef.Sho,(stef.W), '-.k');
            %             plot(tits.Sho,(tits.E), '-b');
            %             plot(tits.Sho,(tits.W), '-.b');
            %             hold off
            %         end
        end
        
        
        %--------------------------------------------------
        % plot for paper with Stefan and Yuri
        %--------------------------------------------------
        if (paper==1)
            
            close all
            
            crop = 1;
            
            lim.xmax = lim.xmax - eps_val * (lim.xmax-lim.xmin) *dt;
            lim.zmax = lim.zmax + eps_val * (lim.zmax-lim.zmin) *dt;
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            szzd  = hdf5read(filename,'/Centers/sxxd'); szzd = cast(szzd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            
            
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T      = reshape(T,params(4)-1,params(5)-1)'- 273;
            
            
            VizGrid = PhaseMap( filename, VizGrid  );
            % Interpolate composition on strainrate grid
            Phase = griddata(VizGrid.x_plot,VizGrid.z_plot', VizGrid.ph,xc_plot,zc_plot');
            
            %         Efact_max = max(max(eII(Phase==0)))/eps_val
            %
            %         A = zeros(size(eII));
            %         A(Phase==0) = eII(Phase==0)/eps_val;
            %         figure(10)
            %         hold on
            %         pcolor(xc_plot, zc_plot, A)
            %         shading interp
            %         colorbar
            %         hold off
            
            
            xc  = xc_plot;
            zc  = zc_plot;
            if crop == 1
                [Phase,  xc_plot, zc_plot] = CropCellArray( Phase,  xc, zc, lim );
                [eII,  xc_plot, zc_plot] = CropCellArray( eII,  xc, zc, lim );
                [T,    ~, ~] = CropCellArray(   T,  xc, zc, lim );
            end
            
            
            if show == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off'), clf
            end
            
            if show == 1
                % eII
                subplot (121)
                k = imagesc(xc_plot,zc_plot,log10(eII));
                colorbar('Location', 'SouthOutside')
                if exist('minE', 'var') caxis([minE maxE]); end
                shading flat; axis xy image
                hold on;
            end
            
            [c, hb] = contour(xc_plot,zc_plot,log10(1/1*eII./eps_val), [0.0 0.0]);
            set(hb, 'Color', 'w', 'LineWidth', 2)
            
            if size(c,2)>0
                
                %             [a,trans] = find(c(2,:)>1.5*lim.zmax)
                [a,trans] = find(c(1,:)==0);
                
                if size(trans,2) == 2
                    
                    % find the lower line
                    cont1.x = c(1,trans(1)+1:trans(2)-1);
                    cont1.z = c(2,trans(1)+1:trans(2)-1);
                    P1 = polyfit(cont1.x, cont1.z, 1);
                    if show == 1; plot(xc_plot, P1(1).*xc_plot + P1(2), '-k'); end;
                    
                    % find the lower line
                    cont2.x = c(1,trans(2)+1:end);
                    cont2.z = c(2,trans(2)+1:end);
                    P2 = polyfit(cont2.x, cont2.z, 1);
                    if show == 1; plot(xc_plot, P2(1).*xc_plot + P2(2), '-k'); end;
                    
                    % lower line point
                    EqLow = P2(1).*xc_plot + P2(2);
                    LowInsideWindow = EqLow < max(zc_plot) & EqLow > min(zc_plot);
                    [u, v] = find( LowInsideWindow );
                    xLow = mean( xc_plot(u) );
                    zLow = mean( EqLow(u) );
                    
                    % upper line point
                    EqUp = P1(1).*xc_plot + P1(2);
                    UpInsideWindow = EqUp < max(zc_plot) & EqUp > min(zc_plot);
                    [u, v] = find( UpInsideWindow );
                    xUp = mean( xc_plot(u) );
                    zUp = mean( EqUp(u) );
                    
                    if show == 1;  plot(xLow, zLow, '-r*'); end;
                    if show == 1;  plot(xUp, zUp, '-r*'); end;
                    
                    
                    % Find the middle line
                    a = mean([P1(1), P2(1)]);
                    b = mean([P1(2), P2(2)]);
                    if show == 1;  plot(xc_plot, a.*xc_plot + b, '-.w'); end;
                    
                    %% Find the normal equation
                    
                    
                    
                    % Point on the center line
                    Pz = 0.5*(zUp + zLow);
                    Px = (Pz-b)/a;
                    
                    
                    
                    if show == 1;  plot(Px, Pz, '-w*'); end;
                    
                    %                 slope = (zUp-zLow) / (xUp-xLow);
                    %                 ord   = zUp-slope*xUp;
                    
                    %                 Normal = slope.*(xSpace) + ord;
                    
                    %                 mean_x = 0.5*(max(xc_plot) + min(xc_plot)); %mean_x = 0.5*(max(xSpace) + min(xSpace))
                    %                 mean_z = 0.5*(min(a.*xc_plot + b) + max(a.*xc_plot + b)); mean_z = 0.5*(min(a.*xSpace + b) + max(a.*xSpace + b));
                    %
                    anorm  = -1/a;
                    bnorm  =  (Px*a + b) - anorm*Px;
                    
                    % Intersect normal and contour
                    Xi1 = -(P2(2) - bnorm) / (P2(1) - anorm);
                    Xi2 = -(P1(2) - bnorm) / (P1(1) - anorm);
                    Xi = min(Xi1, Xi2);
                    Zi = anorm*Xi + bnorm;
                    
                    % create space
                    (pi-atan(a))*90/pi;
                    atan(a)
                    xmin   = Xi;
                    Width  = 1.5 * 1/1 * ( abs(P1(2)-P2(2))  *cos(atan(a)));%sqrt( (xUp - xLow).^2 + (zUp - zLow).^2 );
                    dx     = Width / 101;
                    xmin   = xmin;% - 0.1*xmin;
                    xSpace = xmin:dx:xmin+1/4*(xc_plot(end)-xc_plot(1));%+Width;
                    
                    Normal = -1/a.*(xSpace) + bnorm;
                    
                    if show == 1;  plot(Xi, Zi, '-w*'); end
                    %                 xSpace = xc_plot;
                    
                    if show == 1;  plot(xSpace, Normal, '-r', 'LineWidth', 2); end;
                    %                 xii = xc_plot;
                    xii = xSpace;
                    zii = Normal;
                    dist = sqrt((xii-xii(1)).^2 + (zii-zii(1)).^2 );
                    Eii_itp = griddata(xc_plot,zc_plot,eII,xii,zii);
                    Eii_itp = Eii_itp-min(Eii_itp);
                    Eii_itp = Eii_itp(isnan(Eii_itp)==0);
                    xii = xii(isnan(Eii_itp)==0);
                    zii = zii(isnan(Eii_itp)==0);
                    dist = sqrt((xii-xii(1)).^2 + (zii-zii(1)).^2 );
                    
                    % find max of Eii_itp
                    [~, ie]=max(Eii_itp);
                    Amp = max(Eii_itp);
                    fit = diff(log(Eii_itp)) ./ diff((dist-dist(ie)).^2);
                    
                    % ------------------------
                    
                    if show == 1
                        
                        for i = 1:1
                            if ( sum (sum (VizGrid.ph == i-1)) ~= 0)
                                Vizbuf = VizGrid.ph;
                                Vizbuf(VizGrid.ph == i) = 1;
                                Vizbuf(VizGrid.ph ~= i) = 69;
                                [c2, hb] = contour(VizGrid.x_plot, VizGrid.z_plot, Vizbuf, [ 1 1]);
                                set(hb, 'Color', 'k', 'LineWidth', 2);
                            end
                        end
                        
                        hold off
                        xlabel(xLabel, 'Fontsize', Ftsz), ylabel(zLabel,'Fontsize', Ftsz);
                        %                     title(['eII at' TimeLabel, ' min = ', num2str(min(eII(:)),'%2.2e'), ' max = ', num2str(max(eII(:)))], 'Fontsize', Ftsz)
                        if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
                        
                    end
                    
                    % ------------------------
                    minThick = dx/2;
                    maxThick = 0.25*(xg_coord(end)- xg_coord(1));
                    dThick   = (maxThick-minThick)/500;
                    fit = [minThick:dThick:maxThick];
                    if L > 1e3; fit = fit./1e3; end
                    fit = -fit;
                    
                    % Find best fit
                    error = zeros(size(fit));
                    for o=1:length(error)
                        subplot (122)
                        hold on
                        plot(dist-dist(ie), Eii_itp, dist-dist(ie), Amp*exp(-1/(fit(o))^2*(dist-dist(ie)).^2))
                        %                     model    = Amp*exp(-1/(fit(o))^2 * (dist-dist(ie)).^2);
                        model    = Amp*exp(-1/2/(fit(o))^2 * (dist-dist(ie)).^2);
                        %                     sqrt( sum(abs(model - Eii_itp).^2) )
                        error(o) = sqrt( sum(abs(model - Eii_itp).^2) );
                    end
                    
                    [~, ibest] = min(error);
                    bestfit = fit(ibest);
                    
                    if show == 1
                        subplot(122)
                        hold on
                        plot(dist-dist(ie), Eii_itp, '.k', 'MarkerSize', 10)
                        %                     plot(dist-dist(ie), Amp*exp(-1/(bestfit)^2*(dist-dist(ie)).^2), '-.g', 'LineWidth', 2)
                        plot(dist-dist(ie), Amp*exp(-1/2/(bestfit)^2*(dist-dist(ie)).^2), '-.g', 'LineWidth', 2)
                        
                        axis auto
                        hold off
                        box off
                        drawnow
                        set(gca, 'fontsize', 15);
                        %                     axis([-8 8 0 5e-13])
                        %                     set(gca, 'Ytick', [0 1 2 3 4 5]'.*1e-13, 'YtickLabel', ['0'; '1e-13'; '2e-13' ;'3e-13' ;'4e-13'; '5e-13']);
                        % %                     set(gca, 'Xtick', [-8 -4 0 4 8].*1);
                        %                     set(gca, 'position', [0.1 0.1 0.25 0.25].*[1 21/29.7 1 21/29.7])
                        %                     orient tall
                        
                    end
                    
                    SbThickYuri(icount) = -bestfit;
                    if (xg_coord(end) - xg_coord(1)) > 1e3
                        SbThickYuri(icount) = SbThickYuri(icount)*1e3;
                    end
                    
                    
                    % Slope
                    %                 min_b = min(P1(2), P2(2));
                    %                 if (min_b==P1(2))
                    %                     angle = atand(P1(1));
                    %                 elseif (min_b==P2(2))
                    %                     angle = atand(P2(1));
                    %                 end
                    angle = atand(a);
                    
                    
                    % Shearband angle
                    SbAngle(icount) = angle;
                    
                    % Shearband thickness
                    mean_x = mean(xc_plot);
                    h1 = P1(1)*mean_x + P1(2);
                    h2 = P2(1)*mean_x + P2(2);
                    H  = abs(h1-h2);
                    
                    SbThick(icount) = H*sind(90-SbAngle(icount));
                    if (xg_coord(end) - xg_coord(1)) > 1e3
                        SbThick(icount) = SbThick(icount)*1e3;
                    end
                end
            end
            
            % Temperature increase
            if icount == 1
                Tref =  T(1,1);
            end
            
            %         figure(90)
            %         subplot 121
            %         k = imagesc(xc_plot,zc_plot,T);
            %         axis xy image
            %         colorbar
            %          subplot 122
            %         k = imagesc(xc_plot,zc_plot,log10(eII));
            %         axis xy image
            %         colorbar
            %         if exist('minE', 'var') caxis([minE maxE]); end
            
            Tchar = T(end,1);
            SbDT(icount)  = max(T(:)) - Tref;
            SbDTY(icount) = max(T(:)) - Tchar;
            
            % Strain rate amplification
            SbEiiFact(icount) = max(eII(:))/eps_val;
            
            % Shortening/time
            short(icount)  = strain;
            short2(icount) = strain2;
            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/Fig_StrainRateCompo',num2str(istep,'%05d'),file_suffix], format, res)
                
            end
            
            if istep == iend
                save ('TransientData', 'short', 'SbThick', 'SbEiiFact', 'SbDTY', 'SbDT', 'SbAngle')
            end
            close all
            clf
        end
        
        %--------------------------------------------------
        % plot for paper with Stefan and Yuri
        %--------------------------------------------------
        if (paper2==1)
            
            FTSZ = 18;
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            szzd  = hdf5read(filename,'/Centers/sxxd'); szzd = cast(szzd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T      = reshape(T,params(4)-1,params(5)-1)'- 273;
            
            p_vec_x1 = zeros(size(eII));
            p_vec_y1 = zeros(size(eII));
            p_vec_x2 = zeros(size(eII));
            p_vec_y2 = zeros(size(eII));
            
            angle = deg2rad(45);
            rot = [cos(angle) -sin(angle ); sin(angle) cos(angle)];
            
            
            for i=1:size(eII,1)
                for j=1:size(eII,2)
                    tau = [exx(i,j) exz(i,j); exz(i,j) -exx(i,j)];
                    [V,Q]   = eig(tau);
                    
                    rot     = [cos(angle) -sin(angle ); sin(angle) cos(angle)];
                    Q       = rot * Q * rot'
                    V(:,1)  = rot*V(:,1);
                    V(:,2)  = rot*V(:,2);
                    
                    p_vec_x1(i,j) = V(1,1);
                    p_vec_y1(i,j) = V(2,1);
                    p_vec_x2(i,j) = V(1,2);
                    p_vec_y2(i,j) = V(2,2);
                    
                    
                end
            end
            
            
            
            if crop == 1
                [eII,  xc_plot, zc_plot] = CropCellArray( eII,  xc_plot, zc_plot, lim );
                [T,    xc_plot, zc_plot] = CropCellArray(   T,  xc_plot, zc_plot, lim );
            end
            
            VizGrid = PhaseMap( filename, VizGrid  );
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount),clf
            else
                figure('Visible', 'Off')
            end
            
            colormap('Bone')        % eII
            
            set(gca, 'fontsize', FTSZ)
            
            hold on
            p = imagesc(xc_plot*1e0,zc_plot*1e0,log10(eII));
            %         colorbar('FontSize', 15)
            colorbar('Location', 'SouthOutside', 'FontSize', FTSZ)
            if exist('minE', 'var') caxis([minE maxE]); end
            
            %          [c, hb] = contour(xc_plot*1e0,zc_plot*1e0,log10(eII/eps_val), [0 0]);
            %          set(hb, 'Color', 'w', 'LineWidth', 2 )
            %
            for i = 1:1
                if ( sum (sum (VizGrid.ph == i-1)) ~= 0)
                    Vizbuf = VizGrid.ph;
                    Vizbuf(VizGrid.ph == i) = 1;
                    Vizbuf(VizGrid.ph ~= i) = 69;
                    [c2, hb] = contour(VizGrid.x_plot, VizGrid.z_plot, Vizbuf, [ 1 1]);
                    set(hb, 'Color', 'k', 'LineWidth', 1.5);
                end
            end
            
            hold off
            shading flat
            alpha(p, 0.65)
            xlabel(xLabel, 'FontSize', FTSZ), ylabel(zLabel, 'FontSize', FTSZ);
            title(['eII at' TimeLabel, ' min = ', num2str(min(eII(:)),'%2.2e'), ' max = ', num2str(max(eII(:)))],'FontSize', 15)
            axis xy image
            drawnow
            
            xlim([-35 35]); ylim([-25 25])
            
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/Fig_StrainRateCompoBIG',num2str(istep,'%05d'),file_suffix], format, res)
                close all
                
            end
        end
        
        
        %--------------------------------------------------
        % Kinectic energy
        %--------------------------------------------------
        if ( kinetic_e == 1 )
            
            rho_s = hdf5read(filename,'/Vertices/rho_s');
            rho_s = cast(rho_s, 'double');
            VxNodes.Vx = hdf5read(filename,'/VxNodes/Vx'); VxNodes.Vx = cast(VxNodes.Vx, 'double');
            VzNodes.Vz = hdf5read(filename,'/VzNodes/Vz'); VzNodes.Vz = cast(VzNodes.Vz, 'double');
            rho = reshape(rho_s,nx,nz)';
            Vx  = reshape(VxNodes.Vx,nx,nz+1)';
            Vz  = reshape(VzNodes.Vz,nx+1,nz)';
            
            % Cell centered velocity
            Vx   = Vx(2:end-1,:);
            Vz   = Vz(:,2:end-1);
            Vxc  = 0.5 * (Vx(:,2:end) + Vx(:,1:end-1));
            Vzc  = 0.5 * (Vz(2:end,:) + Vz(1:end-1,:));
            rhoc = 0.25*(rho(1:end-1,1:end-1) + rho(1:end-1, 2:end) + rho(2:end, 1:end-1) + rho(2:end,2:end));
            
            V = (xg_coord(end)-xg_coord(1))*(zg_coord(end)-zg_coord(1));
            
            E_kin(icount) =  1/2 * sum(sum(rhoc.*(Vxc.^2 + Vzc.^2))) *dx*dz; % 1/V
            t_vec(icount) = time;
            
            if istep>0
                figure(90)
                hold on
                plot(t_vec(icount), E_kin(icount), '.')
            end
            
        end
        
        %--------------------------------------------------
        % plot stress and strain rate invariants
        %--------------------------------------------------
        if ( strainR_div == 1 )
            
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            ezz = -exx;
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            
            eII_Vx = 0.5*(eII(:,1:end-1) + eII(:,2:end)); eII_Vx = eII_Vx(2:end-1,:);
            eII_Vz = 0.5*(eII(1:end-1,:) + eII(2:end,:)); eII_Vz = eII_Vz(:,2:end-1);
            
            % dEii/dx + dEii/dz
            divEii = sqrt(1/dx * diff(eII_Vx, 1, 2).^2 + 1/dz * diff(eII_Vz, 1, 1).^2 );
            
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            imagesc( xc_plot(2:end-1), zc_plot(2:end-1), (divEii) )
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['eII at' TimeLabel, ' max = ', num2str((max(eII(:))), '%2.2e'  ), ' s^{-1}' ])
            %         if exist('minE', 'var') caxis([minE maxE]); end
            %         caxis([0 2e-15])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/Fig_Stress_StrainRdiv',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
        end
        
        %--------------------------------------------------
        % plot stress and strain rate invariants
        %--------------------------------------------------
        if ( Christmas_Tree == 1 )
            
            x1 = xc_coord(1) + 0/2*abs(xc_coord(end)-xc_coord(1));
            x2 = xc_coord(1) + 10/10*abs(xc_coord(end)-xc_coord(1));
            
            [aa,b1] = min(abs(xc_coord-x1));
            [aa,b2] = min(abs(xc_coord-x2));
            
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)'-273;
            P = hdf5read(filename,'/Centers/P'); P = cast(P, 'double');
            P = reshape(P,params(4)-1,params(5)-1)';
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            eta = 0.25*(eta(1:end-1,1:end-1) + eta(2:end,1:end-1) + eta(1:end-1,2:end) + eta(2:end,2:end));
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            ezz = -exx;
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            
            sxxd = (reshape(sxxd,params(4)-1,params(5)-1)');
            szzd = -sxxd;
            sxz  = (reshape(sxz, params(4)  ,params(5)  )');
            sII  = sqrt( 0.5*(2*sxxd.^2 + 0.5*(sxz(1:end-1,1:end-1).^2 + sxz(2:end,1:end-1).^2 + sxz(1:end-1,2:end).^2 + sxz(2:end,2:end).^2 ) ) );
            
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            HSc  = 2*sxzc.*exzc + sxxd.*exx + szzd.*ezz;
            %         minHS   = min(HSc(:));
            %         maxHS   = max(HSc(:));
            %         minSxxd = min(sxxd(:))
            %         maxSxxd = max(sxxd(:))
            %         minTxz  = min(sxz(:))
            %         maxTxz  = max(sxz(:))
            %         min(sII(:))
            %         max(sII(:))
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            % sII
            subplot(141)
            hold on
            plot( (sII(:,b1)), zc_plot, '-b'  )
            plot( (sII(:,b2)), zc_plot, '--r'  )
            hold off
            if exist('minSii', 'var') axis([minSii maxSii min(zc_plot) max(zc_plot)]); end
            xlabel(['S_{II} [Mpa]'])
            ylabel(['z [km]'])
            %         legend('1/4', '3/4', 'Location', 'SouthOutside')
            %         xlabel(xLabel), ylabel(zLabel);
            %         title(['sII at' TimeLabel, ' min = ', num2str((min(sII(:))), '%2.2e' ), ' Pa', ' max = ', num2str((max(sII(:))), '%2.2e' ), ' Pa' ])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minSii', 'var') caxis([minSii maxSii]); end
            
            subplot(142)
            hold on
            plot( log10(eta(:,b1)), zc_plot, '-b'  )
            plot( log10(eta(:,b2)), zc_plot, '--r'  )
            hold off
            %         axis([18 25 min(zc_plot) max(zc_plot)])
            xlabel(['\eta [Pa.s]'])
            ylabel(['z [km]'])
            time
            title(['sII at' TimeLabel, ' min = ', num2str((min(sII(:))), '%2.2e' ), ' Pa', ' max = ', num2str((max(sII(:))), '%2.2e' ), ' Pa' ])
            
            %         title(['\eta at' TimeLabel, ' min = ', num2str((min(eta(:))), '%2.2e' ), ' Pa.s', ' max = ', num2str((max(eta(:))), '%2.2e' ), ' Pa.s' ])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(143)
            hold on
            plot( (T(:,b1)), zc_plot, '-b'  )
            plot( (T(:,b2)), zc_plot, '--r'  )
            hold off
            %         axis([0 1500 min(zc_plot) max(zc_plot)])
            xlabel(['T [C]'])
            ylabel(['z [km]'])
            
            subplot(144)
            hold on
            plot( (P(:,b1))/1e6, zc_plot, '-b'  )
            plot( (P(:,b2))/1e6, zc_plot, '--r'  )
            hold off
            %         axis([0 4000 min(zc_plot) max(zc_plot)])
            xlabel(['P [MPa]'])
            ylabel(['z [km]'])
            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/Fig_Profiles',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
            %         save('Yoann_VP_initial', 'zc_plot', 'P', 'T', 'eta', 'sII')
        end
        
        %--------------------------------------------------
        %  Topography
        %--------------------------------------------------
        if ( topo == 1 )
            
            if mdoodz6==1
                x_mark = hdf5read(filename,'/Topo/x_mark');
                z_mark = hdf5read(filename,'/Topo/z_mark');
                Vx_mark = hdf5read(filename,'/Topo/Vx_mark');
                Vz_mark = hdf5read(filename,'/Topo/Vz_mark');
                Vx_grid = hdf5read(filename,'/Topo/Vx_grid');
                Vz_grid = hdf5read(filename,'/Topo/Vz_grid');
                x_mark = cast(x_mark, 'double');
                z_mark = cast(z_mark, 'double');
                height = hdf5read(filename,'/Topo/z_grid');
            else
                height = hdf5read(filename,'/Topo/height');
            end
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            subplot(311)
            title(TimeLabel)
%             title(['min = ', num2str(min(ztopo)), ' max = ', num2str(max(ztopo)), TimeLabel])
            hold on
            plot( x_mark, z_mark-0*mean(z_mark), '.r' )
%             plot( xg_coord, height-0*mean(height), '-o')
            hold off
            
            subplot(312),hold on
            plot( x_mark, Vz_mark, '.r' )
            plot( xvz_coord, Vz_grid, '-o')
            subplot(313),hold on
            plot( x_mark, Vx_mark, '.r' )
            
            
            % figure(45)
            % plot( time/365.25/3600/24/1e6, max(ztopo), '+')
            % axis([0 20 0 900])
            % xlabel('Time [ky]')
            % ylabel('max. Topo [km]')
            % hold on
            
%             mxtopo(icount) = max(ztopo);
            drawnow
            
            if istep == iend
                save( [path,'TopoData'], 'mxtopo', 't_vec')
                disp('file saved')
            end
            
            if printfig == 1
                print([path, './Fig_Topo/Fig_Topo', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % topography and viscosity plot
        %--------------------------------------------------
        if (topo_eta_plot==1)
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta_n = hdf5read(filename,'/Centers/eta_n');  eta_n = cast(eta_n, 'double');
            
            eta   = reshape(eta_s,params(4),params(5))';
            eta_n = reshape(eta_n,params(4)-1,params(5)-1)';
            
            xtopo = hdf5read(filename,'/Topo/x_mark');
            ztopo = hdf5read(filename,'/Topo/z_mark');
            xtopo = cast(xtopo, 'double');
            ztopo = cast(ztopo, 'double');
            height = hdf5read(filename,'/Topo/z_grid');
            
            % Nodes
            if print2screen == 1
                figCount = figCount +1;
                f=figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            subplot(211)
            
            title(['min = ', num2str(min(ztopo)), ' max = ', num2str(max(ztopo)), TimeLabel])
            hold on
            plot( xtopo, ztopo, '.r' )
            plot( xg_coord, height, '-o')
            hold off
            
            subplot(212)
            hold on
            imagesc(xg_plot, zg_plot, log10(eta));
            shading flat, axis xy image, colorbar;
            title(['Viscosity (log10) at' TimeLabel])
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            if Ccontours == 1; imagesc(xg_plot, zg_plot, log10(eta), 'Visible', 'Off'); end
            
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minEta', 'var') caxis([minEta maxEta]); end
            drawnow
            
            if printfig == 1
                print([path, './Fig_Viscosity/Fig_TopoViscosity', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % topography and strain rate
        %--------------------------------------------------
        if (topo_SR_plot==1)
            
            
            xtopo = hdf5read(filename,'/Topo/x_mark');
            ztopo = hdf5read(filename,'/Topo/z_mark');
            xtopo = cast(xtopo, 'double');
            ztopo = cast(ztopo, 'double');
            height = hdf5read(filename,'/Topo/z_grid');
            
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            
            % Nodes
            if print2screen == 1
                figCount = figCount +1;
                f=figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            subplot(211)
            
            title(['min = ', num2str(min(ztopo)), ' max = ', num2str(max(ztopo)), TimeLabel])
            hold on
            plot( xtopo, ztopo, '.r' )
            plot( xg_coord, height, '-o')
            hold off
            
            subplot(212)
            hold on
            imagesc(xg_plot, zg_plot, log10(eII));
            shading flat, axis xy image, colorbar('SouthOutside');
            title(['Strain rate invariant II (log10) at' TimeLabel])
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            if Ccontours == 1; imagesc(xg_plot, zg_plot, log10(eII), 'Visible', 'Off'); end
            
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minEii', 'var') caxis([minEii maxEii]); end
            drawnow
            
            if printfig == 1
                print([path, './Fig_Stress_StrainRate/Fig_TopoSR', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        %  Topography through time
        %--------------------------------------------------
        if ( topo_maps == 1 )
            
            if mdoodz6 == 1
                vx_grid = hdf5read(filename,'/Topo/Vx_grid');
                vz_grid = hdf5read(filename,'/Topo/Vz_grid');
                vx_grid = cast(vx_grid, 'double');
                vz_grid = cast(vz_grid, 'double');
                z_grid  = hdf5read(filename,'/Topo/z_grid');
                z_grid  = cast(z_grid, 'double');
            else
                vx_grid = hdf5read(filename,'/Topo/vx');
                vz_grid = hdf5read(filename,'/Topo/vz');
                vx_grid = cast(vx_grid, 'double');
                vz_grid = cast(vz_grid, 'double');
                z_grid  = hdf5read(filename,'/Topo/height');
                z_grid  = cast(z_grid, 'double');
            end

            if istep == istart
                topo2D = zeros(iend/ijump+1, length(xg_coord));
                vx2D   = zeros(iend/ijump+1, length(xg_coord));
                vz2D   = zeros(iend/ijump+1, length(xg_coord));
                x2D    = zeros(iend/ijump+1, length(xg_coord));
                t2D    = zeros(iend/ijump+1, length(xg_coord));
                t2Dc   = zeros(iend/ijump+1, length(xg_coord)-1);
            end
            
            topo2D(icount,:) = z_grid;
            x2D(icount,:)    = xg_coord;
            vx2D(icount,:)   = vx_grid;
            vz2D(icount,:)   = 0.5*(vz_grid(1:end-1)+vz_grid(2:end-0));
            t2D(icount,:)    = time;
            t2Dc(icount,:)   = time;
            
            if istep == iend
                if print2screen == 1
                    figCount = figCount +1;
                    figure(figCount), clf
                else
                    figure('Visible', 'Off')
                end
                
                subplot(221)
                pcolor(x2D/1e3,t2D/My,topo2D/1e3-topo_ref), colorbar, shading flat;
                caxis([-3 2])
                xlabel('x [km]'), ylabel('t [My]')
                title('Topography [km]')
                if exist('minHkm', 'var'), caxis([minHkm maxHkm]); end

                
                subplot(223)
                pcolor(x2D/1e3,t2D/My,vx2D*y*1e3), colorbar, shading flat;
                
                xlabel('x [km]'), ylabel('t [My]')
                title('Surface velocity x [mm/y]')
                if exist('minVxmmy', 'var'), caxis([minVxmmy maxVxmmy]); end
                
                subplot(224)
                pcolor(x2D/1e3,t2D/My,vz2D*y*1e3), colorbar, shading flat;
                hold on 
                if exist('minVymmy', 'var'), caxis([minVymmy maxVymmy]); end
%                 contour(x2D/1e3,t2D/My,topo2D/1e3-topo_ref, [0 0], 'k')

                xlabel('x [km]'), ylabel('t [My]')
                title('Surface velocity z [mm/y]')
                
                
                dt2D     = t2D(2:end,:) - t2D(1:end-1,:);
                dx2D     = x2D(:,2:end) - x2D(:,1:end-1);
                Vz_wrong = (topo2D(2:end,:) - topo2D(1:end-1,:)) ./ dt2D;
                t2D_av   = 0.5*(t2D(2:end,:) + t2D(1:end-1,:));
                x2D_av   = 0.5*(x2D(:,2:end) + x2D(:,1:end-1));
                vx2D_av   = 0.5*(vx2D(:,2:end) + vx2D(:,1:end-1));
                
                dHdx     = (topo2D(:,2:end) - topo2D(:,1:end-1)) ./ dx2D;
                
                Vz_wrong_av = 0.5*(Vz_wrong(:,2:end) + Vz_wrong(:,1:end-1));
                Vz_corr  = Vz_wrong_av + vx2D_av(1:end-1,:).*dHdx(1:end-1,:);
                
                 subplot(222)
                pcolor(x2D_av(1:end,:)/1e3,t2Dc/My,log10(abs(diff(vx2D,1,2)/(dx)))), colorbar, shading flat;
%                 caxis([-2.5 2.5])
                xlabel('x [km]'), ylabel('t [My]')
                title('Exx [1/s]')
                
                
%                 subplot(235)
%                 pcolor(x2D(1:end-1,:)/1e3,t2D_av/My,Vz_wrong*y*1e3), colorbar, shading flat;
% %                 caxis([-2.5 2.5])
%                 xlabel('x [km]'), ylabel('t [My]')
%                 title('DH/Dt = dH/dt [mm/y]')
%                 
%                 subplot(236)
%                 pcolor(x2D_av(1:end-1,:)/1e3,t2D_av(:,1:end-1)/My,Vz_corr*y*1e3), colorbar, shading flat;
% %                 caxis([-2.5 2.5])
%                 xlabel('x [km]'), ylabel('t [My]')
%                 title('DH/Dt = dH/dt + Vx*dHdx [mm/y]')
                
                if printfig == 1
                    print([path, './Fig_TopoMaps', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end

            end
        end
        
        %--------------------------------------------------
        % phases + uplift rate
        %--------------------------------------------------
        if(phases_uplift==1)
            
            if mdoodz6 == 1
                vx_grid = hdf5read(filename,'/Topo/Vx_grid');
                vz_grid = hdf5read(filename,'/Topo/Vz_grid');
                vx_grid = cast(vx_grid, 'double');
                vz_grid = cast(vz_grid, 'double');
                z_grid  = hdf5read(filename,'/Topo/z_grid');
                z_grid  = cast(z_grid, 'double');
            else
                vx_grid = hdf5read(filename,'/Topo/vx');
                vz_grid = hdf5read(filename,'/Topo/vz');
                vx_grid = cast(vx_grid, 'double');
                vz_grid = cast(vz_grid, 'double');
                z_grid  = hdf5read(filename,'/Topo/height');
                z_grid  = cast(z_grid, 'double');
            end
            
            VizGrid = PhaseMap_hr( filename, VizGrid );
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)';

            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            clf
            
%             subplot(211)
%             hold on
%             plot(xvz_coord/1e3, vz_grid*y*1e3, '-k', 'LineWidth', 2)
%             plot(xvz_coord/1e3, 0*vz_grid*y*1e3, '-.k')
% 
%             
%             subplot(212)
% 
%             hold on
            
%             colormap(map);
%             imagesc(VizGrid.x_plot_hr, VizGrid.z_plot_hr, VizGrid.ph_hr);
%             set(gca,'Ydir','Normal')
%             shading interp, caxis ([-1 size(map,1)-2])%, axis ij
%             axis image
%             title(['Time = ', num2str(time/1e6),' Myrs'], 'FontSize', 15, 'FontName', 'Myriad pro');
             contour(VizGrid.x_plot_hr, VizGrid.z_plot_hr, VizGrid.ph_hr, 1:10, 'k')
            
%             % Isotherms
%             hold on
%             v = 200:200:2000;
%             [c, h] = contour(VizGrid.x_plot, VizGrid.z_plot, T-273, v);
%             set(h, 'Color', 'w', 'LineWidth', 1.0);
%             hold off


            hold on
            vz_grid_av = 0.5*(vz_grid(1:end-1)+vz_grid(2:end));
%             plot(xg_coord(vz_grid_av>0)/1e3, z_grid(vz_grid_av>0)/1e3, '.r', 'LineWidth', 4)
%             plot(xg_coord(vz_grid_av<0)/1e3, z_grid(vz_grid_av<0)/1e3, '.b', 'LineWidth', 4)
            
            scatter(xg_coord/1e3, z_grid/1e3, 60, vz_grid_av*y*1e3, 'filled')
            colormap('jet')
            colorbar
            caxis([-2.5 2.5])
            hold off
            axis image
            
%             colorbar('Location', 'SouthOutside')
            xlabel(xLabel, 'FontSize', Ftsz), ylabel(zLabel, 'FontSize', Ftsz)
            
            title(['Surface uplift at' TimeLabel], 'FontSize', Ftsz);
            % %         hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            xlim([min(VizGrid.x_plot) max(VizGrid.x_plot)]); ylim([-40 10])
            
            %         if exist('minPhase', 'var') caxis([minPhase maxPhase]); end
%             set(gca, 'FontSize', Ftsz);%, 'FontName', 'Myriad pro')
            set(gcf,'Color','white')
            
            if printfig == 1
                if crop == 1
                    print([path, './Fig_Phases/Crop_PhasesUpliftGrid', num2str(istep,'%05d'),file_suffix], format, res)
                else
                    print([path, './Fig_Phases/Fig_PhasesUpliftGrid', num2str(istep,'%05d'),file_suffix], format, res)
                end
                close all
            end
            
            clear VizGrid.ph VizGrid.x VizGrid.z
            
        end
       
        %--------------------------------------------------
        % plot phases on visualisation grid
        %--------------------------------------------------
        if (srate_add_perc==1)
            
            eII_el  = hdf5read(filename,'/Centers/eII_el');
            eII_el  = cast(eII_el , 'double');
            eII_el = reshape(eII_el,params(4)-1,params(5)-1)';
            
            eII_pl  = hdf5read(filename,'/Centers/eII_pl');
            eII_pl  = cast(eII_pl , 'double');
            eII_pl = reshape(eII_pl,params(4)-1,params(5)-1)';
            
            eII_pwl  = hdf5read(filename,'/Centers/eII_pwl');
            eII_pwl  = cast(eII_pwl , 'double');
            eII_pwl = reshape(eII_pwl,params(4)-1,params(5)-1)';
            
            eII_exp  = hdf5read(filename,'/Centers/eII_exp');
            eII_exp  = cast(eII_exp , 'double');
            eII_exp = reshape(eII_exp,params(4)-1,params(5)-1)';
            
            eII_lin  = hdf5read(filename,'/Centers/eII_lin');
            eII_lin  = cast(eII_lin , 'double');
            eII_lin = reshape(eII_lin,params(4)-1,params(5)-1)';
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            colormap('jet');
            
            subplot(321)
            imagesc( xc_plot, zc_plot, eII_el./eII )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (eII_el./eII), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            caxis([0 1]);
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel)
            title(['Elastic strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(323)
            imagesc( xc_plot, zc_plot, eII_pl./eII )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (eII_pl./eII), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            caxis([0 1]);
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Plastic strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(322)
            imagesc( xc_plot, zc_plot, eII_pwl./eII )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, log10(eII_pwl./eII), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            caxis([0 1]);
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Power law strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(324)
            imagesc( xc_plot, zc_plot, eII_exp./eII )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, ((eII_exp./eII)), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            caxis([0 1]);
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Exponential strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(325)
            imagesc( xc_plot, zc_plot, eII_lin./eII )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, ((eII_lin./eII)), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            caxis([0 1]);
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Linear strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(326)
            imagesc( xc_plot, zc_plot, eII./eII )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, ((eII./eII)), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            caxis([0 1]);
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Total strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/StrainRatesAdd',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot phases on visualisation grid
        %--------------------------------------------------
        if(phase_temp==1)
            
            VizGrid = PhaseMap( filename, VizGrid );
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)';

%             if print2screen == 1
%                 figCount = figCount +1;
%                 fh=figure(figCount); clf;
%             else
%                 fh=figure('Visible', 'Off');
%                 set(fh,'units','normalized','position',[0.0 0.0 1.0 1.0]);
%             end
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end

            hold on
           
            clf
            colormap(map);
            imagesc(VizGrid.x_plot, VizGrid.z_plot, VizGrid.ph);
            set(gca,'Ydir','Normal')
            shading interp, caxis ([-1 size(map,1)-2])%, axis ij
            axis image
            title(['Time = ', num2str(time/1e6),' Myrs'], 'FontSize', 15, 'FontName', 'Myriad pro');
            
            
            % Isotherms
            hold on
            v = 200:200:2000;
            [c, h] = contour(VizGrid.x_plot, VizGrid.z_plot, T-273, v);
            set(h, 'Color', 'w', 'LineWidth', 1.0);
            hold off
            
            colorbar('Location', 'SouthOutside')
            xlabel(xLabel, 'FontSize', Ftsz), ylabel(zLabel, 'FontSize', Ftsz)
            title(['Phases at' TimeLabel], 'FontSize', Ftsz);
            % %         hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minPhase', 'var') caxis([minPhase maxPhase]); end
            set(gca, 'FontSize', Ftsz);%, 'FontName', 'Myriad pro')
            set(gcf,'Color','white')
            
            if printfig == 1
                if crop == 1
                    print([path, './Fig_Phases/Crop_PhasesGrid', num2str(istep,'%05d'),file_suffix], format, res)
                else
                    print([path, './Fig_Phases/Fig_PhasesGrid', num2str(istep,'%05d'),file_suffix], format, res)
                end
                close all
            end
            
            if printfig==2
                
                gifname = [path, './Fig_Phases/Crop_PhasesGrid.gif'];
                axis tight manual % this ensures that getframe() returns a consistent size
                % Capture the plot as an image
                frame = getframe(fh);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                % Write to the GIF File
                
                %             [imind(:,:,1,icount),cm]=rgb2ind(im,256);
                if icount == 1
                    imwrite(imind,cm,gifname,'gif', 'Loopcount',inf);
                else
                    imwrite(imind,cm,gifname,'gif','WriteMode','append');
                end
            end

            clear VizGrid.ph VizGrid.x VizGrid.z
            
        end
        
        %--------------------------------------------------
        % plot phases on visualisation grid
        %--------------------------------------------------
        if(phase_temp2==1)
            
            VizGrid = PhaseMap_hr( filename, VizGrid );
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)';

            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end

            hold on
            
            clf
            colormap(map);
            imagesc(VizGrid.x_plot_hr, VizGrid.z_plot_hr, VizGrid.ph_hr);
            set(gca,'Ydir','Normal')
            shading interp, caxis ([-1 size(map,1)-2])%, axis ij
            axis image
            title(['Time = ', num2str(time/1e6),' Myrs'], 'FontSize', 15, 'FontName', 'Myriad pro');
            
            
            % Isotherms
            hold on
            v = 200:200:2000;
            [c, h] = contour(VizGrid.x_plot, VizGrid.z_plot, T-273, v);
            set(h, 'Color', 'w', 'LineWidth', 1.0);
            hold off
            
            colorbar('Location', 'SouthOutside')
            xlabel(xLabel, 'FontSize', Ftsz), ylabel(zLabel, 'FontSize', Ftsz)
            title(['Phases at' TimeLabel], 'FontSize', Ftsz);
            % %         hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minPhase', 'var') caxis([minPhase maxPhase]); end
            set(gca, 'FontSize', Ftsz);%, 'FontName', 'Myriad pro')
            set(gcf,'Color','white')
            
            if printfig == 1
                if crop == 1
                    print([path, './Fig_Phases/Crop_PhasesGrid', num2str(istep,'%05d'),file_suffix], format, res)
                else
                    print([path, './Fig_Phases/Fig_PhasesGrid', num2str(istep,'%05d'),file_suffix], format, res)
                end
                close all
            end
            
            if printfig==2
                
                gifname = [path, './Fig_Phases/Crop_PhasesGrid.gif'];
                axis tight manual % this ensures that getframe() returns a consistent size
                % Capture the plot as an image
                frame = getframe(fh);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                % Write to the GIF File
                
                %             [imind(:,:,1,icount),cm]=rgb2ind(im,256);
                if icount == 1
                    imwrite(imind,cm,gifname,'gif', 'Loopcount',inf);
                else
                    imwrite(imind,cm,gifname,'gif','WriteMode','append');
                end
            end

            clear VizGrid.ph VizGrid.x VizGrid.z
            
        end
        
        %--------------------------------------------------
        % plot grain size
        %--------------------------------------------------
        if(grain_size==1)
            
            T = hdf5read(filename,'/Centers/T');    T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)';
            d = hdf5read(filename,'/Centers/d');    d = cast(d, 'double');
            d = reshape(d,params(4)-1,params(5)-1)';
            eta_n = hdf5read(filename,'/Centers/eta_n'); eta_n = cast(eta_n, 'double');
            eta_n = reshape(eta_n,params(4)-1,params(5)-1)';
            
            eII_el  = hdf5read(filename,'/Centers/eII_el');
            eII_el  = cast(eII_el , 'double');
            eII_el = reshape(eII_el,params(4)-1,params(5)-1)';
            
            eII_pl  = hdf5read(filename,'/Centers/eII_pl');
            eII_pl  = cast(eII_pl , 'double');
            eII_pl = reshape(eII_pl,params(4)-1,params(5)-1)';
            
            eII_pwl  = hdf5read(filename,'/Centers/eII_pwl');
            eII_pwl  = cast(eII_pwl , 'double');
            eII_pwl = reshape(eII_pwl,params(4)-1,params(5)-1)';
            
            eII_exp  = hdf5read(filename,'/Centers/eII_exp');
            eII_exp  = cast(eII_exp , 'double');
            eII_exp = reshape(eII_exp,params(4)-1,params(5)-1)';
            
            eII_lin  = hdf5read(filename,'/Centers/eII_lin');
            eII_lin  = cast(eII_lin , 'double');
            eII_lin = reshape(eII_lin,params(4)-1,params(5)-1)';
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            
            exx  = (reshape(exxd,params(4)-1,params(5)-1)');
            exz  = (reshape(exz, params(4)  ,params(5)  )');
            txx  = (reshape(sxxd,params(4)-1,params(5)-1)');
            txz  = (reshape(sxz, params(4)  ,params(5)  )');
            exzn = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            txzn = 0.25*(txz(1:end-1,1:end-1) + txz(2:end,1:end-1) + txz(1:end-1,2:end) + txz(2:end,2:end));
            %         eII  = sqrt( exx.^2 + exzn.^2);
            %         tII  = sqrt( txx.^2 + txzn.^2);
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            tII = sqrt( 0.5*( 2*txx.^2 + 0.5*(txz(1:end-1,1:end-1).^2 + txz(2:end,1:end-1).^2 + txz(1:end-1,2:end).^2 + txz(2:end,2:end).^2 ) ) );
            Hs  = 4*eta_n.*( exzn.^2 + exx.^2);
            
            % Covey-Crump deformation maps
            pg              = 3.0;
            Kg              = 2.5e9* power(10,-6*pg);
            Qg              = 175e3                 ;
            gam             = 1                     ;
            cg              = pi;
            lambda          = 0.1;
            R               = 8.314;
            
            d1 = 1e-6:5e-7:5e-3;
            T1 = (200:600)+273;
            one_T1 = 1./T1;
            size(d1)
            size(one_T1)
            [one_T2,d2] = meshgrid(one_T1,d1);
            diss = (Kg.*exp(-Qg./R.*one_T2)*cg*gam)./(pg*lambda.*d2.^(1+pg));
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            colormap('jet')
            imagesc(xc_plot, zc_plot, log10(d));
            shading flat% caxis ([-1 size(map,1)-2])%, axis ij
            axis image
            %         colorbar('Location', 'SouthOutside')
            xlabel(xLabel, 'FontName', 'Myriad pro'), ylabel(zLabel, 'FontName', 'Myriad pro')
            caxis([-5 -2.75])
            title(['Grain size at' TimeLabel], 'FontSize', 12, 'FontName', 'Myriad pro');
            set(gca, 'FontSize', 12, 'FontName', 'Myriad pro'), colorbar;%, 'FontName', 'Myriad pro')
            
            %         if print2screen == 1
            %             figCount = figCount +1;
            %             figure(figCount), clf
            %         else
            %             figure('Visible', 'Off')
            %         end
            %
            %         subplot(2,2,[1 3])
            %         hold on
            %
            % %         ic = double(int16(nx/2));
            % %         jc = double(int16(nz/2));
            % %         dpinch = d(jc,ic);
            % %         dswell = d(jc,1);
            % %         plot(350, (dpinch*1e6), 'b+', 'MarkerSize', 15)
            % %         plot(350, (dswell*1e6), 'r+', 'MarkerSize', 15)
            %         ilay    = find(d(:)>1.1e-5);
            %         GSlayer = d(ilay);
            %         dmax    = max(GSlayer);
            %         dmin    = min(GSlayer);
            %         plot(350, (dmin*1e6), 'b+', 'MarkerSize', 15)
            %         plot(350, (dmax*1e6), 'r+', 'MarkerSize', 15)
            %
            % %         [C,h] = contour(one_T2*1000,log(d2),log10(diss*1e-6), [-16:2:-4],'-k');
            % %         [C,h] = contour(one_T2*1000,(d2*1e6),log10(diss), [-10:2:0],'-k');
            %          [C,h] = contour(1./one_T2-273,(d2*1e6),log10(diss), [-10:2:0],'-k');
            %
            %
            %         set(h,'ShowText','on','TextStep',get(h,'LevelStep')*1)
            %         hold off
            %         hl=legend('min. Layer D', 'max. Layer D')
            %         set(hl, 'box', 'off', 'Location', 'NorthWest')
            %         axis([200 500  8 5000])%
            % %         title('log Diss [MPa/s]', 'FontSize', 12, 'FontName', 'Myriad pro')
            %         title('log Diss [Pa/s]', 'FontSize', 12, 'FontName', 'Myriad pro')
            %
            % %         ylabel('ln D [\mu m]','FontName', 'Myriad pro')
            %         ylabel(' D [\mum]','FontName', 'Myriad pro')
            %
            % %         xlabel('1/T*10^{3} [K]','FontName', 'Myriad pro')
            %         xlabel('T [C]','FontName', 'Myriad pro')
            %
            %         set(gca, 'FontSize', 12, 'FontName', 'Myriad pro')
            %         set(gca, 'xscale', 'log')
            %         set(gca, 'yscale', 'log')
            %         set(gca,'Ytick',[1 10 100 1000 5000],'YTickLabel',{'1', '10', '100', '1000', '5000'})
            %         set(gca,'Xtick',[200 300 400 500],'XTickLabel',{'200', '300', '400', '500'})
            %
            %         subplot(2,2,2)
            %         imagesc(xc_plot, zc_plot, log10(Hs));
            %         shading flat% caxis ([-1 size(map,1)-2])%, axis ij
            %         axis image
            % %         title(['Time = ', num2str(time/1e6),' Myrs'], 'FontSize', 15, 'FontName', 'Myriad pro');
            % %         colorbar('Location', 'SouthOutside')
            % %         xlabel(xLabel)
            % %         caxis([-16 -4])
            %         caxis([-10 0])
            %         ylabel(zLabel, 'FontName', 'Myriad pro')
            %         title(['log Diss at' TimeLabel], 'FontSize', 12, 'FontName', 'Myriad pro');
            %         set(gca, 'FontSize', 12, 'FontName', 'Myriad pro'), colorbar;%, 'FontName', 'Myriad pro')
            %
            % %         subplot (312)
            % %         imagesc(xc_plot, zc_plot, eII_pwl./eII);
            % %         shading flat% caxis ([-1 size(map,1)-2])%, axis ij
            % %         axis image
            % % %         colorbar('Location', 'SouthOutside')
            % % %         xlabel(xLabel),
            % %         ylabel(zLabel, 'FontName', 'Myriad pro')
            % %         title(['\beta at' TimeLabel], 'FontSize', 12, 'FontName', 'Myriad pro');
            % %         set(gca, 'FontSize', 12, 'FontName', 'Myriad pro'), colorbar;%, 'FontName', 'Myriad pro')
            %
            %
            %         subplot(2,2,4)
            %          imagesc(xc_plot, zc_plot, log10(d*1e6));
            %         shading flat% caxis ([-1 size(map,1)-2])%, axis ij
            %         axis image
            % %         colorbar('Location', 'SouthOutside')
            %         xlabel(xLabel, 'FontName', 'Myriad pro'), ylabel(zLabel, 'FontName', 'Myriad pro')
            %         caxis([0 3])
            %         title(['Grain size at' TimeLabel], 'FontSize', 12, 'FontName', 'Myriad pro');
            %         set(gca, 'FontSize', 12, 'FontName', 'Myriad pro'), colorbar;%, 'FontName', 'Myriad pro')
            %
            %
            % %         % recalc grain size based on piezometer
            % %         d2 = (cg*gam*Kg*exp(-Qg/R/(350+273)) ./ (pg*lambda*eII_pwl.*tII) ).^(1./(1+pg));
            % %         figure(4)
            % %         imagesc(xc_plot, zc_plot, log10(d2));
            % %         shading flat% caxis ([-1 size(map,1)-2])%, axis ij
            % %         axis image
            % % %         colorbar('Location', 'SouthOutside')
            % %         xlabel(xLabel, 'FontName', 'Myriad pro'), ylabel(zLabel, 'FontName', 'Myriad pro')
            % %         caxis([-6 -3])
            % %         title(['Grain size at' TimeLabel], 'FontSize', 12, 'FontName', 'Myriad pro');
            % %         set(gca, 'FontSize', 12, 'FontName', 'Myriad pro'), colorbar;%, 'FontName', 'Myriad pro')
            
            
            
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if printfig == 1
                if crop == 1
                    print([path, './Fig_GS/Crop_GS', num2str(istep,'%05d'),file_suffix], format, res)
                else
                    print([path, './Fig_GS/Fig_GS', num2str(istep,'%05d'),file_suffix], format, res)
                end
                close all
            end
        end
        
        %--------------------------------------------------
        % plot grain size
        %--------------------------------------------------
        if(grain_size_evol==1)
            
            T = hdf5read(filename,'/Centers/T');    T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)';
            d = hdf5read(filename,'/Centers/d');    d = cast(d, 'double');
            d = reshape(d,params(4)-1,params(5)-1)';
            eta_n = hdf5read(filename,'/Centers/eta_n'); eta_n = cast(eta_n, 'double');
            eta_n = reshape(eta_n,params(4)-1,params(5)-1)';
            
            eII_el  = hdf5read(filename,'/Centers/eII_el');
            eII_el  = cast(eII_el , 'double');
            eII_el = reshape(eII_el,params(4)-1,params(5)-1)';
            
            eII_pl  = hdf5read(filename,'/Centers/eII_pl');
            eII_pl  = cast(eII_pl , 'double');
            eII_pl = reshape(eII_pl,params(4)-1,params(5)-1)';
            
            eII_pwl  = hdf5read(filename,'/Centers/eII_pwl');
            eII_pwl  = cast(eII_pwl , 'double');
            eII_pwl = reshape(eII_pwl,params(4)-1,params(5)-1)';
            
            eII_exp  = hdf5read(filename,'/Centers/eII_exp');
            eII_exp  = cast(eII_exp , 'double');
            eII_exp = reshape(eII_exp,params(4)-1,params(5)-1)';
            
            eII_lin  = hdf5read(filename,'/Centers/eII_lin');
            eII_lin  = cast(eII_lin , 'double');
            eII_lin = reshape(eII_lin,params(4)-1,params(5)-1)';
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            
            exx  = (reshape(exxd,params(4)-1,params(5)-1)');
            exz  = (reshape(exz, params(4)  ,params(5)  )');
            txx  = (reshape(sxxd,params(4)-1,params(5)-1)');
            txz  = (reshape(sxz, params(4)  ,params(5)  )');
            exzn = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            txzn = 0.25*(txz(1:end-1,1:end-1) + txz(2:end,1:end-1) + txz(1:end-1,2:end) + txz(2:end,2:end));
            %         eII  = sqrt( exx.^2 + exzn.^2);
            %         tII  = sqrt( txx.^2 + txzn.^2);
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            tII = sqrt( 0.5*( 2*txx.^2 + 0.5*(txz(1:end-1,1:end-1).^2 + txz(2:end,1:end-1).^2 + txz(1:end-1,2:end).^2 + txz(2:end,2:end).^2 ) ) );
            Hs  = 4*eta_n.*( exzn.^2 + exx.^2);
            
            eff = eII_pwl./eII; eff(eff>1.0) = 1.0;
            
            % Where to pick data
            ic = double(int16(nx/2));
            jc = double(int16(nz/2));
            dpinch(icount) = d(jc,ic)*1e6;
            dswell(icount) = d(jc,1)*1e6;
            
            tIIpinch(icount) = tII(jc,ic)*1e-6;
            tIIswell(icount) = tII(jc,1)*1e-6;
            
            Edispinch(icount) = eff(jc,ic);
            Edisswell(icount) = eff(jc,1);
            
            tnd(icount)  = params(1)*1e-12;
            
            %          if print2screen == 1
            %             figCount = figCount +1;
            %             figure(figCount), %clf
            %         else
            %             figure('Visible', 'Off')
            %          end
            %
            %           subplot(1,3,1)
            %          imagesc(xc_plot, zc_plot, tII);
            %         shading flat% caxis ([-1 size(map,1)-2])%, axis ij
            %         axis image
            % %         colorbar('Location', 'SouthOutside')
            %         xlabel(xLabel, 'FontName', 'Myriad pro'), ylabel(zLabel, 'FontName', 'Myriad pro')
            % %         caxis([17 22])
            %         title(['Viscosity  at' TimeLabel], 'FontSize', 12, 'FontName', 'Myriad pro');
            %         set(gca, 'FontSize', 12, 'FontName', 'Myriad pro'), colorbar;%, 'FontName', 'Myriad pro')
            %
            %          subplot(1,3,2)
            %          imagesc(xc_plot, zc_plot, eff);
            %         shading flat% caxis ([-1 size(map,1)-2])%, axis ij
            %         axis image
            % %         colorbar('Location', 'SouthOutside')
            %         xlabel(xLabel, 'FontName', 'Myriad pro'), ylabel(zLabel, 'FontName', 'Myriad pro')
            % %         caxis([0 1])
            %         title(['Grain size  at' TimeLabel], 'FontSize', 12, 'FontName', 'Myriad pro');
            %         set(gca, 'FontSize', 12, 'FontName', 'Myriad pro'), colorbar;%, 'FontName', 'Myriad pro')
            %
            %
            % %          subplot(1,3,3)
            %          imagesc(xc_plot, zc_plot, log10(d));
            %         shading flat% caxis ([-1 size(map,1)-2])%, axis ij
            %         axis image
            % %         colorbar('Location', 'SouthOutside')
            %         xlabel(xLabel, 'FontName', 'Myriad pro'), ylabel(zLabel, 'FontName', 'Myriad pro')
            %         caxis([-5 -2.8])
            %         title(['Grain size  at' TimeLabel], 'FontSize', 12, 'FontName', 'Myriad pro');
            %         set(gca, 'FontSize', 12, 'FontName', 'Myriad pro'), colorbar;%, 'FontName', 'Myriad pro')
            
            %         % recalc grain size based on piezometer
            %         % Covey-Crump deformation maps
            %         pg              = 3.0;
            %         Kg              = 2.5e9* power(10,-6*pg);
            %         Qg              = 175e3                 ;
            %         gam             = 1                     ;
            %         cg              = pi;
            %         lambda          = 0.1;
            %         R               = 8.314;
            %         d2 = (cg*gam*Kg*exp(-Qg/R/(350+273)) ./ (pg*lambda*eII_pwl.*tII) ).^(1./(1+pg));
            %         d2(d<2e-5) = 1e-5;
            %         figure(4)
            %         imagesc(xc_plot, zc_plot, log10(d2));
            %         shading flat% caxis ([-1 size(map,1)-2])%, axis ij
            %         axis image
            % %         colorbar('Location', 'SouthOutside')
            %         xlabel(xLabel, 'FontName', 'Myriad pro'), ylabel(zLabel, 'FontName', 'Myriad pro')
            %         caxis([-5 -2.8])
            %         title(['Piezometer Grain size check at' TimeLabel], 'FontSize', 12, 'FontName', 'Myriad pro');
            %         set(gca, 'FontSize', 12, 'FontName', 'Myriad pro'), colorbar;%, 'FontName', 'Myriad pro')
            
            
            
            %           if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if printfig == 1
            %             if crop == 1
            %                 print([path, './Fig_GS/Crop_GS', num2str(istep,'%05d'),file_suffix], format, res)
            %             else
            %                 print([path, './Fig_GS/Fig_GS', num2str(istep,'%05d'),file_suffix], format, res)
            %             end
            %             close all
            %         end
            
            %         figure(43)
            %         subplot(131), hold on
            %         plot(tnd(icount), tIIpinch(icount), 'r+', tnd(icount),  tIIswell(icount), 'b+')
            %         xlabel( 't^*')
            %         ylabel( '\tau_{II} [MPa]')
            %
            %         subplot(132), hold on
            %         plot(tnd(icount),   dpinch(icount), 'r+', tnd(icount),    dswell(icount), 'b+')
            %         xlabel( 't^*')
            %         ylabel( 'D [\mum]')
            %
            %         subplot(133), hold on
            %         plot(tnd(icount), Edispinch(icount), 'r+', tnd(icount), Edisswell(icount), 'b+')
            %         xlabel( 't^*')
            %         ylabel( 'Edisl/Etot')
            
            %         subplot(131), hold on
            %         plot(strain, tIIpinch(icount), 'rx', strain,  tIIswell(icount), 'b+')
            %         xlabel( '\epsilon %')
            %         ylabel( '\tau_{II} [MPa]')
            %
            %         subplot(132), hold on
            %         plot(strain,   dpinch(icount), 'rx', strain,    dswell(icount), 'b+')
            %         xlabel( '\epsilon %')
            %         ylabel( 'D [\mum]')
            %
            %         subplot(133), hold on
            %         plot(strain, Edispinch(icount), 'rx', strain, Edisswell(icount), 'b+')
            %         xlabel( '\epsilon %')
            %         ylabel( '\beta')
            
        end
        
        if(gs_stefan==1)
            
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)'-273;
            d     = hdf5read(filename,'/Centers/d');    d = cast(d, 'double');
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            eta = 0.25*(eta(1:end-1,1:end-1) + eta(2:end,1:end-1) + eta(1:end-1,2:end) + eta(2:end,2:end));
            
            d   = (reshape(d   ,params(4)-1,params(5)-1)');
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            ezz = -exx;
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            %         sII = 2*eII.*eta;
            
            sxxd = (reshape(sxxd,params(4)-1,params(5)-1)');
            szzd = -sxxd;
            sxz  = (reshape(sxz, params(4)  ,params(5)  )');
            sII  = sqrt( 0.5*(2*sxxd.^2 + 0.5*(sxz(1:end-1,1:end-1).^2 + sxz(2:end,1:end-1).^2 + sxz(1:end-1,2:end).^2 + sxz(2:end,2:end).^2 ) ) );
            
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            HSc = 2*sxzc.*exzc + sxxd.*exx + szzd.*ezz;
            
            VizGrid = PhaseMap( filename, VizGrid );
            
            if crop == 1
                [ VizGrid.ph, VizGrid.x_plot, VizGrid.z_plot ] = CropCellArray( VizGrid.ph, VizGrid.x_plot, VizGrid.z_plot, lim );
            end
            
            lay_thick            = sum(VizGrid.ph,1)*dz;
            neck_stress(icount)  = sII(fix(nz/2),fix(nx/2));
            swell_stress(icount) = sII(fix(nz/2),1);
            neck_gs(icount)      = d(fix(nz/2),fix(nx/2));
            swell_gs(icount)     = d(fix(nz/2),1);
            neck_srate(icount)   = eII(fix(nz/2),fix(nx/2));
            swell_srate(icount)  = eII(fix(nz/2),1);
            neck_thick(icount)   = min(lay_thick);
            swell_thick(icount)  = max(lay_thick);
            neck_T(icount)       = T(fix(nz/2),fix(nx/2));
            swell_T(icount)      = T(fix(nz/2),1);
            strainvec(icount)    = strain;
            
            %         if print2screen == 1
            %             figCount = figCount +1;
            %             figure(figCount), clf
            %         else
            %             figure('Visible', 'Off')
            %         end
            %          imagesc(xc_plot, zc_plot, HSc);
            %         shading flat, axis xy image, colorbar;
            %         title(['x  centers at' TimeLabel])
            %         xlabel(xLabel), ylabel(zLabel);
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minX', 'var') caxis([minX maxX]); end
            
            %         plot(strainvec(1:icount),log10(neck_srate(1:icount)) , '.')
            % %         colormap(map);
            % %         p = imagesc(VizGrid.x_plot, VizGrid.z_plot, VizGrid.ph);
            % %         xlabel(xLabel), ylabel(zLabel)
            % %         title(['Phases at' TimeLabel]);
            % %         shading flat, axis xy image
            % %         hold off
            % %         colorbar
            % %         colormap('gray')
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minPhase', 'var') caxis([minPhase maxPhase]); end
            %
            %
            %
            %         grid(gca, 'minor')
            %         alpha (p,0.7)
            %
            %         if printfig == 1
            %             if crop == 1
            %                 print([path, './Fig_Phases/Crop_PhasesGrid', num2str(istep,'%05d'),file_suffix], format, res)
            %             else
            %                 print([path, './Fig_Phases/Fig_PhasesGrid', num2str(istep,'%05d'),file_suffix], format, res)
            %             end
            %             close all
            %         end
            
            clear VizGrid.ph VizGrid.x VizGrid.z
            
        end
        
        %--------------------------------------------------
        % x plot
        %--------------------------------------------------
        if (X_plot==1)
            x = hdf5read(filename,'/Centers/X');  x = cast(x, 'double');
            
            % Nodes
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            imagesc(xc_plot, zc_plot, (reshape(x,nx-1,nz-1)'));
            shading flat, axis xy image, colorbar;
            title(['x  centers at' TimeLabel])
            xlabel(xLabel), ylabel(zLabel);
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minX', 'var') caxis([minX maxX]); end
            drawnow
            
            if printfig == 1
                print([path, './Fig_X/Fig_X', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot strain rate components
        %--------------------------------------------------
        if ( srate_ratio == 1 )
            
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            ezzd  = -exxd;
            %         ezzd  = hdf5read(filename,'/Centers/ezzd'); ezzd = cast(ezzd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            
            exxd = (reshape(exxd,params(4)-1,params(5)-1)');
            ezzd = (reshape(ezzd,params(4)-1,params(5)-1)');
            exz  = (reshape(exz, params(4)  ,params(5)  )');
            exzn = 0.25*(exz(1:end-1,1:end-1) + exz(1:end-1,2:end) + exz(2:end,1:end-1) + exz(2:end,2:end));
            eiin = sqrt(exxd.^2 + exzn.^2);
            gam  = 1*exzn;
            
            viz = log10(abs(exxd)./abs(exzn));
            viz = zeros(size(exzn));
            indPS = eiin>1e-17 & abs(exxd)>abs(gam);
            indSS = eiin>1e-17 & abs(exxd)<abs(gam);
            ratPS = abs(exxd( indPS ))./abs(gam( indPS ));
            ratSS = abs(gam( indSS ))./abs(exxd( indSS ));
            %         viz( indPS ) = (1 * ratPS ./max(ratPS(:)));
            %         viz( indSS ) =-(1 * ratSS ./max(ratSS(:)));
            viz( indPS ) = (1 );
            viz( indSS ) =-(1 );
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            % Create 13-level blue/red diverging color map:
            level = 13; n = ceil(level/2);
            cmap1 = [linspace(1, 1, n); linspace(0, 1, n); linspace(0, 1, n)]';
            cmap2 = [linspace(1, 0, n); linspace(1, 0, n); linspace(1, 1, n)]';
            cmap = [cmap1; cmap2(2:end, :)];
            
            imagesc(xc_plot, zc_plot, viz );
            shading flat, axis xy image, colorbar;
            title(['Exx / Exz  centers at' TimeLabel])
            xlabel(xLabel), ylabel(zLabel);
            colormap(cmap);
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minX', 'var') caxis([minX maxX]); end
            %          caxis([-0.000001 0.000001]);
            caxis([-1 1]);
            
            
            
            if printfig == 1
                print([path,'./Fig_StrainRateComponents/Fig_srate_ratio',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot strain rate components
        %--------------------------------------------------
        if ( PTmax == 1 )
            
            x1 = 20e3;
            x2 = 20e3;
            
            [aa,b1] = min(abs(xc_coord-x1));
            [aa,b2] = min(abs(xc_coord-x2));
            
            Pmax  = hdf5read(filename,'/Centers/Pmax'); Pmax = cast(Pmax, 'double');
            Tmax  = hdf5read(filename,'/Centers/Tmax'); Tmax = cast(Tmax, 'double');
            x0    = hdf5read(filename,'/Centers/xi'); x0   = cast(x0  , 'double');
            z0    = hdf5read(filename,'/Centers/zi'); z0   = cast(z0  , 'double');
            
            P  = hdf5read(filename,'/Centers/P'); Pmax = cast(P, 'double');
            
            P    = (reshape(P   ,params(4)-1,params(5)-1)');
            Pmax = (reshape(Pmax,params(4)-1,params(5)-1)');
            Tmax = (reshape(Tmax,params(4)-1,params(5)-1)');
            x0   = (reshape(x0  ,params(4)-1,params(5)-1)');
            z0   = (reshape(z0  ,params(4)-1,params(5)-1)');
            
            %         if print2screen == 1
            %             figCount = figCount +1;
            %             figure(figCount), clf
            %         else
            %             figure('Visible', 'Off')
            %         end
            %
            %         subplot (211)
            %         imagesc(xc_plot, zc_plot, Pmax );
            %         shading flat, axis xy image, colorbar;
            %         title(['Pmax  centers at' TimeLabel])
            %         xlabel(xLabel), ylabel(zLabel);
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minP', 'var') caxis([minP maxP]); end
            %
            %         subplot (212)
            %         imagesc(xc_plot, zc_plot,Tmax );
            %         shading flat, axis xy image, colorbar;
            %         title(['Tmax  centers at' TimeLabel])
            %         xlabel(xLabel), ylabel(zLabel);
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minX', 'var') caxis([minX maxX]); end
            %
            %         if print2screen == 1
            %             figCount = figCount +1;
            %             figure(figCount), clf
            %         else
            %             figure('Visible', 'Off')
            %         end
            %
            %         subplot (211)
            %         imagesc(xc_plot, zc_plot, x0 );
            %         shading flat, axis xy image, colorbar;
            %         title(['x0  centers at' TimeLabel])
            %         xlabel(xLabel), ylabel(zLabel);
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minX', 'var') caxis([minX maxX]); end
            %
            %         subplot (212)
            %         imagesc(xc_plot, zc_plot, z0 );
            %         shading flat, axis xy image, colorbar;
            %         title(['z0  centers at' TimeLabel])
            %         xlabel(xLabel), ylabel(zLabel);
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minX', 'var') caxis([minX maxX]); end
            %
            %         if print2screen == 1
            %             figCount = figCount +1;
            %             figure(figCount), clf
            %         else
            %             figure('Visible', 'Off')
            %         end
            %
            %         dPdz = diff(Pmax/1e5/1e3,1,1)/(dz*1e-3);
            %         imagesc(xc_plot, zg_plot(2:end-1), dPdz );
            %         shading flat, axis xy image, colorbar;
            %         title(['dPdz  centers at' TimeLabel])
            %         xlabel(xLabel), ylabel(zLabel);
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         caxis([-1.71 1.71]);
            %
            %         if print2screen == 1
            %             figCount = figCount +1;
            %             figure(figCount), clf
            %         else
            %             figure('Visible', 'Off')
            %         end
            %
            %         subplot (211)
            %         [x1,y1]=meshgrid(xc_plot,zc_plot);
            %         xdiff = (x1-x0/1000);
            %         dxmin=(xdiff(1,1));
            %         dxmax=(xdiff(1,end));
            %         ddx  = (dxmax-dxmin)/(length(xc_plot)-1);
            %         dxcor=(dxmin:ddx:dxmax)';
            %         size(dxcor)
            %         size(xc_plot)
            %         [x2,y2]=meshgrid(dxcor,zc_plot);
            %         imagesc(xc_plot, zc_plot, (x1-x0/1000)-x2 );
            %         shading flat, axis xy image, colorbar;
            %         title(['x0  centers at' TimeLabel])
            %         xlabel(xLabel), ylabel(zLabel);
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %
            %         subplot (212)
            %         [x1,y1]=meshgrid(xc_plot,zc_plot);
            %         imagesc(xg_plot(2:end-1), zc_plot, diff(x1-x0,1,2) );
            %         shading flat, axis xy image, colorbar;
            %         title(['x0  centers at' TimeLabel])
            %         xlabel(xLabel), ylabel(zLabel);
            %         hold on
            %         if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            %         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         caxis([-400 -200])
            
            VizGrid = PhaseMap( filename, VizGrid );
            
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)';
            
            %         if crop == 1
            %             [ VizGrid.ph, VizGrid.x_plot, VizGrid.z_plot ] = CropCellArray( VizGrid.ph, VizGrid.x_plot, VizGrid.z_plot, lim );
            %             [ T, VizGrid.x_plot, VizGrid.z_plot ] = CropCellArray( T, VizGrid.x_plot, VizGrid.z_plot, lim );
            %         end
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            hold on
            
            clf
            colormap(map);
            imagesc(VizGrid.x_plot, VizGrid.z_plot, VizGrid.ph);
            set(gca, 'Ydir', 'normal')
            shading interp, caxis ([-1 size(map,1)-2]), axis image%, axis ij
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            
            figure(56)
            subplot(121)
            %         plot(Pmax(:,b2)/1e9,zc_plot, '-b')
            plot(z0(:,b2)/1e3*-0.29,zc_plot, '-b')
            ylim([lim.zmin lim.zmax])
            subplot(122)
            colormap(map);
            imagesc(VizGrid.x_plot, VizGrid.z_plot, VizGrid.ph);
            set(gca, 'Ydir', 'normal')
            shading interp, caxis ([-1 size(map,1)-2]), axis image%, axis ij
            if crop == 1, xlim([x1/1e3-1 x1/1e3+1]); ylim([lim.zmin lim.zmax]); end
            
            %         plot(VizGrid.ph(:,b2)/1e9,zc_plot, '-b')
            %         ylim([lim.zmin lim.zmax])
            
            
            
            
            if printfig == 1
                print([path,'./Fig_Pressure/Fig_Pmax',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot stress and strain rate invariants
        %--------------------------------------------------
        if ( TminusT0 == 1 )
            
            x1 = xc_coord(1) + 1/2*abs(xc_coord(end)-xc_coord(1));
            x2 = xc_coord(1) + 9/10*abs(xc_coord(end)-xc_coord(1));
            
            [aa,b1] = min(abs(xc_coord-x1));
            [aa,b2] = min(abs(xc_coord-x2));
            
            VizGrid = PhaseMap( filename, VizGrid );
            
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)'-273;
            z0    = hdf5read(filename,'/Centers/zi'); z0   = cast(z0  , 'double');
            z0   = (reshape(z0  ,params(4)-1,params(5)-1)');
            
            if istep == 0
                
                if print2screen == 1
                    figCount = figCount +1;
                    figure(figCount), clf
                else
                    figure('Visible', 'Off')
                end
                
                % Create polynomial fit of temperature profile
                T2 = T(:,b1);
                T3 = T2(T2>eps);
                z3 = zc_plot(T2>eps);
                p  = polyfit(z3,T3, 10);
                hold on
                plot( (T(:,b1)), zc_plot, '-b'  )
                plot( (T(:,b2)), zc_plot, '--r'  )
                plot( polyval(p, z3), z3, '.k'  )
                hold off
                %         axis([0 1500 min(zc_plot) max(zc_plot)])
                xlabel(['T [C]'])
                ylabel(['z [km]'])
            end
            
            % Reconstruction of initial temperature using z0
            T0 = zeros(size(T));
            for icol=1:size(T0,2)
                T0(:,icol) =  polyval( p, z0(:,icol)./1e3);
            end
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            imagesc(xc_plot, zc_plot, T-T0), shading flat, colorbar, axis image
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            hold off
            set(gca, 'Ydir', 'normal')
            title('T-T_0')
            caxis([-200 200])
            
            
            if printfig == 1
                print([path,'./Fig_Temperature/Fig_dT',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
            
            %         save('Yoann_VP_initial', 'zc_plot', 'P', 'T', 'eta', 'sII')
        end
        
        if shear_heating==1
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            eta = reshape(eta_s,params(4),params(5))';
            eta = 0.25*(eta(1:end-1,1:end-1) + eta(2:end,1:end-1) + eta(1:end-1,2:end) + eta(2:end,2:end));
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            ezz = -exx;
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            %         sII = 2*eII.*eta;
            
            sxxd = (reshape(sxxd,params(4)-1,params(5)-1)');
            szzd = -sxxd;
            sxz  = (reshape(sxz, params(4)  ,params(5)  )');
            sII  = sqrt( 0.5*(2*sxxd.^2 + 0.5*(sxz(1:end-1,1:end-1).^2 + sxz(2:end,1:end-1).^2 + sxz(1:end-1,2:end).^2 + sxz(2:end,2:end).^2 ) ) );
            
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            HSc = 2*sxzc.*exzc + sxxd.*exx + szzd.*ezz;
            eII_el  = hdf5read(filename,'/Centers/eII_el');
            eII_el  = cast(eII_el , 'double');
            eII_el = reshape(eII_el,params(4)-1,params(5)-1)';
            HSc(eII_el>0.8*eII) = 1e-18;
            minHS=min(HSc(:))
            maxHS=max(HSc(:))
%             
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            %         if ( H>L ) subplot(121); else subplot(211); end
            
%             p=imagesc( xc_plot, zc_plot, Vmag*(365.25*3600*24*100) );
            p=imagesc( xc_plot, zc_plot, log10(abs(HSc)) );
%             p=imagesc( xc_plot, zc_plot, log10(eII) );

%             alpha(p,0.7);
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end

            imagesc( xc_plot, zc_plot, log10(abs(HSc)) , 'Visible', 'off' );

            shading flat,axis xy image, colorbar% colorbar('SouthOutside');
            
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end

            caxis([-9 -4])
            hold off


            xlabel(xLabel, 'FontSize', 15), ylabel(zLabel, 'FontSize', 15)
            title(['Dissipation at' TimeLabel], 'FontSize', 15)
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            set(gca, 'FontSize', 15);
%            
            
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            
            if printfig == 1
                print([path, './Fig_Temperature/Fig_ShearHeating', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end 
        end
        
        if princi_stress==1
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            P     = hdf5read(filename,'/Centers/P');    P  = cast(P , 'double');
            P = reshape(P,params(4)-1,params(5)-1)';
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            ezz = -exx;
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            %         sII = 2*eII.*eta;
            
            sxxd = (reshape(sxxd,params(4)-1,params(5)-1)');
            szzd = -sxxd;
            sxz  = (reshape(sxz, params(4)  ,params(5)  )');
            sII  = sqrt( 0.5*(2*sxxd.^2 + 0.5*(sxz(1:end-1,1:end-1).^2 + sxz(2:end,1:end-1).^2 + sxz(1:end-1,2:end).^2 + sxz(2:end,2:end).^2 ) ) );
            
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            HSc = 2*sxzc.*exzc + sxxd.*exx + szzd.*ezz;
            eII_el  = hdf5read(filename,'/Centers/eII_el');
            eII_el  = cast(eII_el , 'double');
            eII_el = reshape(eII_el,params(4)-1,params(5)-1)';
            HSc(eII_el>0.8*eII) = 1e-18;
            minHS=min(HSc(:))
            maxHS=max(HSc(:))
%             
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end

            p=imagesc( xc_plot, zc_plot, log10(sII) );
%             alpha(p,0.7);
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end

            imagesc( xc_plot, zc_plot, log10(sII) , 'Visible', 'off' );

            shading flat,axis xy image,  colorbar('SouthOutside');
                        caxis([log10(minSii) log10(maxSii)])

            
            
            % Pricipal stress vectors
            sxx = -P + sxxd;
            szz = -P + szzd;
            p_vec_x1 = zeros(size(eII(1:step:end,1:step:end)));
            p_vec_x3 = zeros(size(eII(1:step:end,1:step:end)));
            p_vec_y1 = zeros(size(eII(1:step:end,1:step:end)));
            p_vec_y3 = zeros(size(eII(1:step:end,1:step:end)));
            s1  = zeros(size(eII(1:step:end,1:step:end)));
            s3  = zeros(size(eII(1:step:end,1:step:end)));
            s1h = zeros(size(eII(1:step:end,1:step:end)));
            
            ii = 0;
            for i=1:step:size(eII,1)
                ii = ii +1;
                jj = 0;
                for j=1:step:size(eII,2)
                    jj  = jj + 1;
                    tau = [sxx(i,j) sxzc(i,j); sxzc(i,j) szz(i,j)];
                    [V,Q]   = eig(tau);
                    p_vec_x3(ii,jj) = V(1,2);
                    p_vec_y3(ii,jj) = V(2,2);
                    p_vec_x1(ii,jj) = V(1,1);
                    p_vec_y1(ii,jj) = V(2,1);
                    s3 (ii,jj) = Q(1,1);
                    s1 (ii,jj) = Q(2,2);
                    s1h(ii,jj) = atand(V(2,1)/V(1,1));
                end
            end
            
            xcq = xc_plot(1:step:end);
            zcq = zc_plot(1:step:end);
%             q=quiver(xcq, zcq,   p_vec_x1,  p_vec_y1, 0.25,'k','Autoscale','on',  'LineWidth', 2, 'ShowArrowHead','off');
            q=quiver(xcq, zcq,  -p_vec_x1, -p_vec_y1, 'w','ShowArrowHead','off',  'LineWidth', 2)
%             q=quiver(xcq, zcq,  -p_vec_x1, -p_vec_y1, 0.25,'k','Autoscale','on', 'LineWidth', 2, 'ShowArrowHead','off');
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end

            hold off

            xlabel(xLabel, 'FontSize', 15), ylabel(zLabel, 'FontSize', 15)
            title(['Principal stresses at' TimeLabel], 'FontSize', 15)
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            set(gca, 'FontSize', 15);
            
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end

            if printfig == 1
                print([path, './Fig_Stress_StrainRate/Fig_PrincipalStress', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end 
        end
            
        
        %--------------------------------------------------
        % plot additive effective strain rates
        %--------------------------------------------------
        if (emergency_benoit==1)
            
            eII_pl  = hdf5read(filename,'/Centers/eII_pl');
            eII_pl  = cast(eII_pl , 'double');
            eII_pl = reshape(eII_pl,params(4)-1,params(5)-1)';
            
            eII_pwl  = hdf5read(filename,'/Centers/eII_pwl');
            eII_pwl  = cast(eII_pwl , 'double');
            eII_pwl = reshape(eII_pwl,params(4)-1,params(5)-1)';
            
            sxxd  = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            sxz   = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz, 'double');
            exxd  = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            exz   = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz, 'double');
            eta_s = hdf5read(filename,'/Vertices/eta_s'); eta_s = cast(eta_s, 'double');
            
            exx = (reshape(exxd,params(4)-1,params(5)-1)');
            exz = (reshape(exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            subplot(311)
            imagesc( xc_plot, zc_plot, (log10(eII)) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, (log10(eII)), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var'); caxis([minEii maxEii]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
                        
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Total strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(312)
            imagesc( xc_plot, zc_plot, log10(eII_pl) )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, log10(eII_pl), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var'); caxis([minEii maxEii]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            caxis([minEii maxEii])
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Frictional strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(313)
            imagesc( xc_plot, zc_plot, log10(eII_pwl) )
            
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, log10(eII_pwl), 'Visible', 'off' )
            shading flat,axis xy image, colorbar;
            if exist('minEii', 'var'); caxis([minEii maxEii]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Ductile strain rate at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            % anti anti-aliasing
            set(gcf,'GraphicsSmoothing','off')
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/StrainRateDucFric',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
        % plot finite strain
        %--------------------------------------------------
        if ( fstrain == 1 )
            
            Fxx  = hdf5read(filename,'/Centers/Fxx'); Fxx = cast(Fxx, 'double');
            Fxx = (reshape(Fxx,params(4)-1,params(5)-1)');
            Fxz  = hdf5read(filename,'/Centers/Fxz'); Fxz = cast(Fxz, 'double');
            Fxz = (reshape(Fxz,params(4)-1,params(5)-1)');
            Fzx  = hdf5read(filename,'/Centers/Fzx'); Fzx = cast(Fzx, 'double');
            Fzx = (reshape(Fzx,params(4)-1,params(5)-1)');
            Fzz  = hdf5read(filename,'/Centers/Fzz'); Fzz = cast(Fzz, 'double');
            Fzz = (reshape(Fzz,params(4)-1,params(5)-1)');
            
            CGxx = zeros(size(Fxx));
            CGxz = zeros(size(Fxx));
            CGzx = zeros(size(Fxx));
            CGzz = zeros(size(Fxx));
            
            for i=1:size(Fxx,2)
                for j=1:size(Fxx,1)
                    
                    F     = [Fxx(j,i) Fxz(j,i); Fzx(j,i) Fzz(j,i)];
                    CG    = F'*F;
                    CGeig = eigs(CG);
                    
                    a=1
                    
                end
            end
            
        end
        
        
        
        
        %--------------------------------------------------
        % plot director vectors
        %--------------------------------------------------
        if(director_vector==1)
            
            VizGrid = PhaseMap_hr( filename, VizGrid );
            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)';
            
            
            ani_fac = hdf5read(filename,'/Centers/ani_fac'); ani_fac = cast(ani_fac, 'double');
            ani_fac = reshape(ani_fac,params(4)-1,params(5)-1)';
            
            
            ndx = hdf5read(filename,'/Centers/nx'); ndx = cast(ndx, 'double'); ndx = reshape(ndx,params(4)-1,params(5)-1)';
            ndz = hdf5read(filename,'/Centers/nz'); ndz = cast(ndz, 'double'); ndz = reshape(ndz,params(4)-1,params(5)-1)';

            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end

            hold on
            
            clf
%             colormap(map);
%             imagesc(VizGrid.x_plot_hr, VizGrid.z_plot_hr, VizGrid.ph_hr);
%             shading interp, caxis ([-1 size(map,1)-2])%, axis ij
            imagesc(xc_plot, zc_plot, ani_fac);
            shading interp,% caxis ([-1 size(map,1)-2])%, axis ij
            set(gca,'Ydir','Normal')
            axis image
            title(['Time = ', num2str(time/1e6),' Myrs'], 'FontSize', 15, 'FontName', 'Myriad pro');
            
           
            % Isotherms
            hold on
            v = 200:200:2000;
            [c, h] = contour(VizGrid.x_plot, VizGrid.z_plot, T-273, v);
            set(h, 'Color', 'w', 'LineWidth', 1.0);
          
            
            % Director
            hv=quiver(xc_plot(1:step:end), zc_plot(1:step:end), ndx(1:step:end,1:step:end), ndz(1:step:end,1:step:end));
%                         hv=quiver(xc_plot(100), zc_plot(100), ndx(100,100), ndz(100,100));
            set(hv, 'Color', 'k', 'LineWidth', 0.5);
            set(hv, 'MaxHeadSize',2.0, 'AutoScaleFactor', 0.5);
           
            hv=quiver(xc_plot(1:step:end), zc_plot(1:step:end), ndz(1:step:end,1:step:end), -ndx(1:step:end,1:step:end));
            set(hv, 'Color', 'w', 'LineWidth', 2);
            set(hv, 'MaxHeadSize',2.0, 'AutoScaleFactor', 0.5, 'ShowArrowHead', 'off');
            hv=quiver(xc_plot(1:step:end), zc_plot(1:step:end), -ndz(1:step:end,1:step:end), ndx(1:step:end,1:step:end));
            set(hv, 'Color', 'w', 'LineWidth', 2);
            set(hv, 'MaxHeadSize',2.0, 'AutoScaleFactor', 0.5, 'ShowArrowHead', 'off');
              hold off
            
            colorbar('Location', 'SouthOutside')
            xlabel(xLabel, 'FontSize', Ftsz), ylabel(zLabel, 'FontSize', Ftsz)
            title(['Phases at' TimeLabel], 'FontSize', Ftsz);
            % %         hold off
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            %         if exist('minPhase', 'var') caxis([minPhase maxPhase]); end
            set(gca, 'FontSize', Ftsz);%, 'FontName', 'Myriad pro')
            set(gcf,'Color','white')
            
            if printfig == 1
                if crop == 1
                    print([path, './Fig_Phases/Crop_PhasesGrid', num2str(istep,'%05d'),file_suffix], format, res)
                else
                    print([path, './Fig_Phases/Fig_PhasesGrid', num2str(istep,'%05d'),file_suffix], format, res)
                end
                close all
            end
            
            if printfig==2
                
                gifname = [path, './Fig_Phases/Crop_PhasesGrid.gif'];
                axis tight manual % this ensures that getframe() returns a consistent size
                % Capture the plot as an image
                frame = getframe(fh);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                % Write to the GIF File
                
                %             [imind(:,:,1,icount),cm]=rgb2ind(im,256);
                if icount == 1
                    imwrite(imind,cm,gifname,'gif', 'Loopcount',inf);
                else
                    imwrite(imind,cm,gifname,'gif','WriteMode','append');
                end
            end

            clear VizGrid.ph VizGrid.x VizGrid.z
            
        end
        
        %--------------------------------------------------
        % plot plastic moduli
        %--------------------------------------------------
        if (Pl_soft==1)
            
            Phi = hdf5read(filename,'/Centers/friction'); Phi = cast(Phi, 'double');
            Phi      = reshape(Phi,params(4)-1,params(5)-1)';
            
            Coh = hdf5read(filename,'/Centers/cohesion'); Coh = cast(Coh, 'double');
            Coh      = reshape(Coh,params(4)-1,params(5)-1)';
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            subplot(211)
            Phi = Phi*180/pi;
            imagesc( xc_plot, zc_plot, Phi )
            shading flat,axis xy image, colorbar;
            colormap('jet')
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, Phi, 'Visible', 'off' )
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Phi at' TimeLabel ' min = ' num2str(min(Phi(:))) ' deg max = ' num2str(max(Phi(:))) ' deg'])
            caxis([15 30])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            subplot(212)
            imagesc( xc_plot, zc_plot, Coh )
            shading flat,axis xy image, colorbar;
            colormap('jet')
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, Coh, 'Visible', 'off' )
            hold off
            xlabel(xLabel), ylabel(zLabel)
            title(['Coh at' TimeLabel ' min = ' num2str(min(Coh(:))) ' MPa max = ' num2str(max(Coh(:))) ' MPa'])
%             caxis([25e6 50e6])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            
            if printfig == 1
                print([path, './Fig_Phi/Fig_Phi', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        
        
        %--------------------------------------------------
        % plot Sole
        %--------------------------------------------------
        if (Sole==1)
            
            %%%%%%%% PHASE
            VizGrid = PhaseMap_hr( filename, VizGrid );

            T = hdf5read(filename,'/Centers/T'); T = cast(T, 'double');
            T = reshape(T,params(4)-1,params(5)-1)' - 273.15;
            
            P  = hdf5read(filename,'/Centers/P');
            P  = cast(P , 'double');
            P = reshape(P,params(4)-1,params(5)-1)';
            
            rho_n = hdf5read(filename,'/Centers/rho_n');
            rho_n = cast(rho_n, 'double');
            rhoc = reshape(rho_n,params(4)-1,params(5)-1)';
            
            % Build lithostatic pressure
            dz = abs(zc_coord(2)-zc_coord(1));
            Plith = P;
            Plith(ncz, :) = rhoc(ncz,:)*gz*dz;
            for o=ncz-1:-1:1
                Plith(o, :) = Plith(o+1, :) - rhoc(o,:)*gz*dz;
            end
            Pdyn = P-Plith;
            
%             if print2screen == 1
%                 figCount = figCount +1;
%                 figure(figCount), clf
%             else
%                 figure('Visible', 'Off')
%             end
%             
%             imagesc( xc_plot, zc_plot, Pdyn )
%             shading flat,axis xy image, colorbar;
%             hold on
%             if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
%             imagesc( xc_plot, zc_plot, Pdyn, 'Visible', 'off' )
%             hold off
%             xlabel(xLabel), ylabel(zLabel)
%             title(['Dynamic pressure at' TimeLabel])
%             if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%             if exist('minPdyn', 'var') caxis([minPdyn maxPdyn]); end
            
            SoleCheck = zeros(size(P));
            SoleCheck( P>0.6e9 & P<1.25e9 & T>700 & T<900) = 1;
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end

        
            clf
            colormap(map);
            imagesc(VizGrid.x_plot_hr, VizGrid.z_plot_hr, VizGrid.ph_hr);
            set(gca,'Ydir','Normal')
            shading interp, caxis ([-1 size(map,1)-2])%, axis ij
            axis image
            title(['Sole ' TimeLabel], 'FontSize', 15, 'FontName', 'Myriad pro');
            hold on    
            [c, h] = contour(VizGrid.x_plot, VizGrid.z_plot, SoleCheck, [0 1]);
            set(h, 'Color', 'w', 'LineWidth', 1.0);
            hold off
            
%             if print2screen == 1
%                 figCount = figCount +1;
%                 figure(figCount), clf
%             else
%                 figure('Visible', 'Off')
%             end
%             
%             imagesc( xc_plot, zc_plot, SoleCheck )
%             shading flat,axis xy image, colorbar;
%             colormap('jet')
%             hold on
%             if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
%             %         imagesc( xc_plot, zc_plot, T, 'Visible', 'off' )
%             
%             imagesc( xc_plot, zc_plot, T, 'Visible', 'off' )
%             hold off
%             
%             xlabel(xLabel), ylabel(zLabel)
%             title(['T at' TimeLabel ' min = ' num2str(min(T(:))) ' C max = ' num2str(max(T(:))) ' C'])
%             if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%             %         if exist('minT', 'var') caxis([minT maxT]); end
%             %         caxis([273 1000])

            if printfig == 1
                print([path, './Fig_Sole', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end
        
        %--------------------------------------------------
    % Overstress plot
    %--------------------------------------------------
    if (overstress==1)
        
        OverS = hdf5read(filename,'/Centers/OverS');
        OverS = cast(OverS, 'double');
        
        min(OverS(:))
        max(OverS(:))
        
        OverS = reshape(OverS,nx-1,nz-1)';
        
        
        load('roma.mat')
        
        if print2screen == 1
            figCount = figCount +1;
            figure(figCount), clf
        else
            figure('Visible', 'Off')
        end
        
        % overstress
        colormap(flipud(roma));
        imagesc( xc_plot, zc_plot, log10(OverS) )
        hold on
        if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
        imagesc( xc_plot, zc_plot, log10(OverS), 'visible', 'off' )
        hold off
        shading flat,axis xy image, colorbar;
        xlabel(xLabel), ylabel(zLabel);
        title(['Overstress at' TimeLabel, ' max = ', num2str((max(OverS(:))), '%2.2e'  ), ' Pa' ])
        if exist('minEii', 'var') caxis([5 log10(3e7)]); end
        if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
        if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
        axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
        
        if printfig == 1
            print([path, './Fig_OverS/Fig_OverS', num2str(istep,'%05d'),file_suffix], format, res)
            close all
        end
        
    end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        if dt_time == 1
            icount
            if icount == 1
                time0 = 0;
            else
                time0 = t_vec(icount-1);
            end
            dt_vec(icount) = (time-time0)/step;
            t_vec(icount)  = time;
            
        end
        
        %     if printfig==2
        %        imwrite(imind,cm,gifname,'gif','Loopcount',inf,'DelayTime', FrameRate)
        %     end
        
        toc
        drawnow
    end
end

if  gs_stefan == 1
    save('GS_transient_1e-14_MR','icount','strainvec','swell_thick','neck_thick','swell_srate','neck_srate','swell_gs','neck_gs','swell_stress','neck_stress','swell_T','neck_T')
end

end

function VizGrid = PhaseMap( filename, VizGrid )

VizGrid.ph = hdf5read(filename,'/VizGrid/compo'); VizGrid.ph = cast(VizGrid.ph, 'double');
VizGrid.ph = (reshape(VizGrid.ph,length(VizGrid.x),length(VizGrid.z))');

% % Clean visualisation grid
% for u = 1:size(VizGrid.ph,1)
%     for v = 1:size(VizGrid.ph,2)
%         if (VizGrid.ph(u,v)==-1)
%
%             % First level cleaning
%             if u>1 && VizGrid.ph(u-1,v)~=-1
%                 VizGrid.ph(u,v) = VizGrid.ph(u-1,v);
%             elseif u<size(VizGrid.ph,1) && VizGrid.ph(u+1,v)~=-1
%                 VizGrid.ph(u,v) = VizGrid.ph(u+1,v);
%             elseif v<size(VizGrid.ph,2) && VizGrid.ph(u,v+1)~=-1
%                 VizGrid.ph(u,v) = VizGrid.ph(u,v+1);
%             elseif v>1 && VizGrid.ph(u,v-1)~=-1
%                 VizGrid.ph(u,v) = VizGrid.ph(u,v-1);
%             end
%
%             if VizGrid.ph(u,v) == -69
%
%                 % Corner level cleaning
%                 if u>1 && v>1 && VizGrid.ph(u-1,v-1)~=-69
%                     VizGrid.ph(u,v) = VizGrid.ph(u-1,v-1);
%                 elseif u<size(VizGrid.ph,1) && v<size(VizGrid.ph,2) && VizGrid.ph(u+1,v+1)~=-69
%                     VizGrid.ph(u,v) = VizGrid.ph(u+1,v+1);
%                 elseif v<size(VizGrid.ph,2) && u>1 && VizGrid.ph(u-1,v+1)~=-69
%                     VizGrid.ph(u,v) = VizGrid.ph(u,v+1);
%                 elseif v>1 && u<size(VizGrid.ph,1) && VizGrid.ph(u+1,v-1)~=-69
%                     VizGrid.ph(u,v) = VizGrid.ph(u,v-1);
%                 end
%             end
%
%             if VizGrid.ph(u,v) == -69
%                 disp('Pretty bad visualition: too few particles in cells !')
%                 % Second level cleaning
%                 if u>2 && VizGrid.ph(u-2,v)~=-69
%                     VizGrid.ph(u,v) = VizGrid.ph(u-2,v);
%                 elseif u<size(VizGrid.ph,1)-1 && VizGrid.ph(u+2,v)~=-69
%                     VizGrid.ph(u,v) = VizGrid.ph(u+2,v);
%                 elseif v<size(VizGrid.ph,2)-1 && VizGrid.ph(u,v+2)~=-69
%                     VizGrid.ph(u,v) = VizGrid.ph(u,v+2);
%                 elseif v>2 && VizGrid.ph(u,v-2)~=-69
%                     VizGrid.ph(u,v) = VizGrid.ph(u,v-2);
%                 end
%             end
%         end
%     end
% end
%
% VizGrid.ph(VizGrid.ph==-69) = 0;
VizGrid.ph(VizGrid.ph==-69) = -1;


end

function VizGrid = PhaseMap_hr( filename, VizGrid )

VizGrid.ph_hr = hdf5read(filename,'/VizGrid/compo_hr'); VizGrid.ph_hr = cast(VizGrid.ph_hr, 'double');
VizGrid.ph_hr = (reshape(VizGrid.ph_hr,length(VizGrid.x_hr),length(VizGrid.z_hr))');

VizGrid.ph_hr(VizGrid.ph_hr==-69) = -1;


end



function PatchCells( ncx, ncy, field, x, y )

imagesc( 0.5*(x(1:end-1) + x(2:end)), 0.5*(y(1:end-1) + y(2:end)), field )

x_tab = zeros(ncx*ncy, 4);
y_tab = zeros(ncx*ncy, 4);
f     = zeros(ncx*ncy,1);
ci = 0;

% Plot cells
for ic = 1:ncx
    for jc = 1:ncy
        ci = ci + 1;
        
        x_arr = [x(ic) x(ic) x(ic+1) x(ic+1)];
        y_arr = [y(jc) y(jc+1) y(jc) y(jc+1)];
        
        x_tab(ci,:) = x_arr;
        y_tab(ci,:) = y_arr;
        
        f(ci) = field(jc,ic);
        
    end
end
id = [ 1 2 4 3];

% patch(x_tab(:,id)', y_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none' );
patch(x_tab(:,id)', y_tab(:,id)', repmat(f,1,4)' );
end


function [ Arr, xc, zc] = CropCellArray( Arr, xc, zc, crop )
[val,icmin] = min(abs(xc-crop.xmin));
[val,icmax] = min(abs(xc-crop.xmax));
[val,jcmin] = min(abs(zc-crop.zmin));
[val,jcmax] = min(abs(zc-crop.zmax));
xc          = xc(icmin:icmax);
zc          = zc(jcmin:jcmax);
Arr         = Arr(jcmin:jcmax,icmin:icmax);
end


function AddCompoContours( filename, VizGrid, crop, lim  )

VizGrid = PhaseMap( filename, VizGrid );

if crop == 1
    [ VizGrid.ph, VizGrid.x_plot, VizGrid.z_plot ] = CropCellArray( VizGrid.ph, VizGrid.x_plot, VizGrid.z_plot, lim );
end

% for i = 0:1:1
%     if ( sum (sum (VizGrid.ph == i-1)) ~= 0)
%         Vizbuf = VizGrid.ph;
%         Vizbuf(VizGrid.ph == i-1) = 1;
%         Vizbuf(VizGrid.ph ~= i-1) = 69;
%         [c, hb] = contour(VizGrid.x_plot, VizGrid.z_plot, Vizbuf, [eps 1+eps]);
%         set(hb, 'Color', 'w', 'LineWidth', 0.5);
%     end
% end

% phases = [-1; 3; 4; 6]
phases = [-1; 0; 3; 4; 6];
for j = 1:length(phases)%max(VizGrid.ph(:))-1
    
        i = phases(j);
        Vizbuf = VizGrid.ph;
        Vizbuf(VizGrid.ph == i) = 1;
        Vizbuf(VizGrid.ph ~= i) = 69;
        [c, hb] = contour(VizGrid.x_plot, VizGrid.z_plot, Vizbuf, [1]);
        set(hb, 'Color', 'w', 'LineWidth', 0.5);
    
end
end