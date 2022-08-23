
function Main

clear
clear all
% close all
clc

DEBUG = 0;

MarkSize=1e0 ;

path = '/Users/romankulakov/CLionProjects/MDOODZ7/cmake-exec/SimpleShearAnisoHomo/'

cd(path)

% File
istart = 0
ijump  = 10;
iend   = 50;

%--------------------------------------------------
% what do you want to plot:
%--------------------------------------------------
stress_evol     = 1;
director_vector = 1;
fstrain         = 1;

data_Stefan     = 0;

if data_Stefan==1
    Stefan = load('../MATLAB_examples/StefanFiniteStrainAnisotropyTest')
    Stefan.Tauxy = [Stefan.Tauxy(1) Stefan.Tauxy];
    Stefan.Nx    = [Stefan.Nx(1) Stefan.Nx];
    Stefan.Ny    = [Stefan.Ny(1) Stefan.Ny];
end


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
single_vector = 1;
step          = 9;
MaskAir       = 1;
ColorFabio    = 1;
y             = 365.25*24*3600;
My            = 1e6*y;
topo_ref      = 0; %6370
mdoodz6       = 1;


% Contour
eps_val = 5e-14;
% eps_val = 1e4;

minRho = 3000;
maxRho = 3600;

% Colorbar limits
minTxx = -11e0;
maxTxx =  0.642;

% minPdyn = -1e9;
% maxPdyn =  1e9;

% minEta = 19;
% maxEta = 25;

% minEii = -16;
% maxEii = -13.5;

minSii = 1e6;
maxSii = 400e6;

% Size of the window
crop       = 1;

lim.xmin   = -0.5;
lim.xmax   =  0.5;
lim.zmin   = -0.5;
lim.zmax   =  0.5;

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
txxvec = zeros(10000,1);
tzzvec = zeros(10000,1);
txzvec = zeros(10000,1);
timevec= zeros(10000,1);
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
    
end

%% LOOP ON OUTPUT FILES
for istep=istart:ijump:iend
    
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
        x_arr = [x_air(:) x_air(:) x_air(:)+dx/1e3 x_air(:)+dx/1e3];
        z_arr = [z_air(:) z_air(:)+dz/1e3 z_air(:) z_air(:)+dz/1e3];
        x_tab(:,:) = x_arr;
        z_tab(:,:) = z_arr;

        %--------------------------------------------------
        % plot stress and strain rate invariants evolution
        %--------------------------------------------------
        if ( stress_evol == 1 )
            
            ndx = hdf5read(filename,'/Centers/nx'); ndx = cast(ndx, 'double'); ndx = reshape(ndx,params(4)-1,params(5)-1)';
            ndz = hdf5read(filename,'/Centers/nz'); ndz = cast(ndz, 'double'); ndz = reshape(ndz,params(4)-1,params(5)-1)';
            
            Cent.sxxd  = hdf5read(filename,'/Centers/sxxd'); Cent.sxxd = cast(Cent.sxxd, 'double');
%             Cent.szzd  = hdf5read(filename,'/Centers/szzd'); Cent.szzd = cast(Cent.szzd, 'double');
            Vert.sxz   = hdf5read(filename,'/Vertices/sxz'); Vert.sxz  = cast(Vert.sxz, 'double');
            Cent.exxd  = hdf5read(filename,'/Centers/exxd'); Cent.exxd = cast(Cent.exxd, 'double');
            Vert.exz   = hdf5read(filename,'/Vertices/exz'); Vert.exz  = cast(Vert.exz, 'double');
            Vert.eta_s = hdf5read(filename,'/Vertices/eta_s'); Vert.eta_s = cast(Vert.eta_s, 'double');
            eta = reshape(Vert.eta_s,params(4),params(5))';
            eta = 0.25*(eta(1:end-1,1:end-1) + eta(2:end,1:end-1) + eta(1:end-1,2:end) + eta(2:end,2:end));
            
            exx = (reshape(Cent.exxd,params(4)-1,params(5)-1)');
            ezz = -exx;
            exz = (reshape(Vert.exz, params(4)  ,params(5)  )');
            eII = sqrt( 0.5*( 2*exx.^2 + 0.5*(exz(1:end-1,1:end-1).^2 + exz(2:end,1:end-1).^2 + exz(1:end-1,2:end).^2 + exz(2:end,2:end).^2 ) ) );
            %         sII = 2*eII.*eta;
            
            sxxd = (reshape(Cent.sxxd,params(4)-1,params(5)-1)');
            szzd = -sxxd;%(reshape(Cent.szzd,params(4)-1,params(5)-1)');
            sxz  = (reshape(Vert.sxz, params(4)  ,params(5)  )');
            sII  = sqrt( 0.5*(2*sxxd.^2 + 0.5*(sxz(1:end-1,1:end-1).^2 + sxz(2:end,1:end-1).^2 + sxz(1:end-1,2:end).^2 + sxz(2:end,2:end).^2 ) ) );
            
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            HSc  = 2*sxzc.*exzc + sxxd.*exx + szzd.*ezz;
            
            sII  = sII(sII>1e-8);
            vol  = length(sII)*dx*dz;
            stress = sum(sII(:).*(dx*dz))/vol;
            if istep == 0
                stress = 0;
            end
            %%%%%%%%%%%%%%%%%%%%%
            figure(90),
            if (istep==0), clf; end
            %%% Stress
            subplot(131), hold on
%             plot(time, mean(sxxd(:)), 'b.')
            plot(time, mean(sxzc(:)), 'ok')
            if istep==iend
                xlabel('time'), title('Txy')
                set(gca, 'FontSize', 16)
            end
            if data_Stefan == 1
                plot(time, Stefan.Tauxy(istep+1), '*r')
            end
            %%% Director
            subplot(132), hold on
            plot(time, mean(ndx(:)), 'ok')
            if istep==iend
                xlabel('time'), title('Nx')
                set(gca, 'FontSize', 16)
            end
            if data_Stefan == 1
                plot(time, Stefan.Nx(istep+1), '*r')
            end
            subplot(133), hold on
            plot(time, mean(ndz(:)), 'ok')
            if data_Stefan == 1
                plot(time, Stefan.Ny(istep+1), '*r')
            end
             if istep==iend
                xlabel('time'), title('Nz')
                set(gca, 'FontSize', 16)
                l=legend('MD6.0', 'Stefan 0D');
                set(l, 'box','off');
            end

            
            eiimaxvec(icount) = max(eII(:));
            short(icount)  =  strain;
            siivec(icount) = sum(sII(:).*(dx*dz))/vol;
            eiivec(icount) = sum(eII(:).*(dx*dz))/vol;
            txxvec(icount) = mean(sxxd(:));
            tzzvec(icount) = mean(szzd(:));
            txzvec(icount) = mean(sxzc(:));
            timevec(icount)= time;
            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/Fig_Stress_StrainRate',num2str(istep,'%05d'),file_suffix], format, res)
                close all
                
            end
            
            if istep == iend
                txxvec(icount+1:end) = [];
                tzzvec(icount+1:end) = [];
                txzvec(icount+1:end) = [];
                timevec(icount+1:end)= [];
%                 save('MaxStressStrainRate', 'siivec', 'eiivec', 'eiimaxvec', 'short')
                 save('StressEvol', 'txxvec', 'tzzvec', 'txzvec', 'timevec')   
            end
            
        end
        
        %--------------------------------------------------
        % plot director vectors
        %--------------------------------------------------
        if (director_vector==1)
            
            VizGrid = PhaseMap_hr( filename, VizGrid );
            Cent.T = hdf5read(filename,'/Centers/T'); Cent.T = cast(Cent.T, 'double');
            T = reshape(Cent.T,params(4)-1,params(5)-1)';
            
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
          
            
            % Director
            if single_vector == 1,
                hv=quiver(xc_plot(1:step:end)*0, zc_plot(1:step:end)*0, ndx(1:step:end,1:step:end), ndz(1:step:end,1:step:end));
            else
                hv=quiver(xc_plot(1:step:end), zc_plot(1:step:end), ndx(1:step:end,1:step:end), ndz(1:step:end,1:step:end));
            end
            set(hv, 'Color', 'k', 'LineWidth', 2.0);
            set(hv, 'MaxHeadSize',2.0, 'AutoScaleFactor', 2/4);
           
            if single_vector == 1,
                hv=quiver(xc_plot(1:step:end)*0, zc_plot(1:step:end)*0, ndz(1:step:end,1:step:end), -ndx(1:step:end,1:step:end));
            else
                hv=quiver(xc_plot(1:step:end), zc_plot(1:step:end), ndz(1:step:end,1:step:end), -ndx(1:step:end,1:step:end));
            end
            set(hv, 'Color', 'w', 'LineWidth', 2.0);
            set(hv, 'MaxHeadSize',2.0, 'AutoScaleFactor', 2/4, 'ShowArrowHead', 'off');
            
              
            % Finite strain
            
            %--------------------------------------------------
            % plot finite strain
            %--------------------------------------------------
            if ( fstrain == 1 )
                
                Cent.Fxx  = hdf5read(filename,'/Centers/Fxx'); Cent.Fxx = cast(Cent.Fxx, 'double');
                Fxx = (reshape(Cent.Fxx,params(4)-1,params(5)-1)');
                Cent.Fxz  = hdf5read(filename,'/Centers/Fxz'); Cent.Fxz = cast(Cent.Fxz, 'double');
                Fxz = (reshape(Cent.Fxz,params(4)-1,params(5)-1)');
                Cent.Fzx  = hdf5read(filename,'/Centers/Fzx'); Cent.Fzx = cast(Cent.Fzx, 'double');
                Fzx = (reshape(Cent.Fzx,params(4)-1,params(5)-1)');
                Cent.Fzz  = hdf5read(filename,'/Centers/Fzz'); Cent.Fzz = cast(Cent.Fzz, 'double');
                Fzz = (reshape(Cent.Fzz,params(4)-1,params(5)-1)');
                
                % Cauchy-Green
                CGxx = Fxx.*Fxx + Fzx.*Fzx;
                CGxz = Fxx.*Fxz + Fzx.*Fzz;
                CGzx = Fxx.*Fxz + Fzx.*Fzz;
                CGzz = Fxz.*Fxz + Fzz.*Fzz;
                
                % Stretch tensor = sqrt(CG)
                [ U0xx,U0xz,U0zx,U0zz]  = Sqrt_2x2( CGxx, CGxz, CGzx, CGzz );
                [e1,e2,Uxx,Uxz,Uzx,Uzz] = Eigs_2x2( U0xx, U0xz, U0zx, U0zz );
                [wxx,wxz,wzx,wzz]       =  Inv_2x2( U0xx, U0xz, U0zx, U0zz );
                [rxx,rxz,rzx,rzz]       =  AxB_2x2( Fxx, Fxz, Fzx, Fzz, wxx, wxz, wzx, wzz );
                [Uxx,Uxz,Uzx,Uzz]       =  AxB_2x2( rxx, rxz, rzx, rzz, Uxx, Uxz, Uzx, Uzz );
                n                       = sqrt(Uxx.^2 + Uzx.^2);
                v11                     = Uxx./n * e1;
                v12                     = Uzx./n * e1;
                n                       = sqrt(Uxz.^2 + Uzz.^2);
                v21                     = Uxz./n * e2;
                v22                     = Uzz./n * e2;
                
                FS_AR = mean(e1(:)) /  mean(e2(:))
                
                % Plot Thibault ellipse
                if single_vector == 1
                    hv = quiver( 0, 0,      mean(v11(:)),      mean(v12(:)) );
                    set(hv, 'Color', 'g', 'LineWidth', 2.0, 'AutoScaleFactor', 1/200, 'ShowArrowHead','off');
                    hv = quiver( 0, 0,     -mean(v11(:)),     -mean(v12(:)) );%,  0.25, 'g', 'AutoScaleFactor', 1/8,  'LineWidth', 3, 'ShowArrowHead','off')
                    set(hv, 'Color', 'g', 'LineWidth', 2.0, 'AutoScaleFactor', 1/200, 'ShowArrowHead','off');
                    hv = quiver( 0, 0,      mean(v21(:)),      mean(v22(:)) );%,  0.25, 'g', 'AutoScaleFactor', 1/8,  'LineWidth', 3, 'ShowArrowHead','off')
                    set(hv, 'Color', 'g', 'LineWidth', 2.0, 'AutoScaleFactor', 1/200, 'ShowArrowHead','off');
                    hv = quiver( 0, 0,     -mean(v21(:)),     -mean(v22(:)) );%,  0.25, 'g', 'AutoScaleFactor', 1/8,  'LineWidth', 3, 'ShowArrowHead','off')
                    set(hv, 'Color', 'g', 'LineWidth', 2.0, 'AutoScaleFactor', 1/200, 'ShowArrowHead','off');
                else
                    hv = quiver( xc_plot(1:step:end), zc_plot(1:step:end),      v11(1:step:end,1:step:end),     v12(1:step:end,1:step:end) );
                    set(hv, 'Color', 'g', 'LineWidth', 2.0, 'AutoScaleFactor', 1/4, 'ShowArrowHead','off');
                    hv = quiver( xc_plot(1:step:end), zc_plot(1:step:end),     -v11(1:step:end,1:step:end),     -v12(1:step:end,1:step:end) );%,  0.25, 'g', 'AutoScaleFactor', 1/8,  'LineWidth', 3, 'ShowArrowHead','off')
                    set(hv, 'Color', 'g', 'LineWidth', 2.0, 'AutoScaleFactor', 1/4, 'ShowArrowHead','off');
                    hv = quiver( xc_plot(1:step:end), zc_plot(1:step:end),      v21(1:step:end,1:step:end),      v22(1:step:end,1:step:end) );%,  0.25, 'g', 'AutoScaleFactor', 1/8,  'LineWidth', 3, 'ShowArrowHead','off')
                    set(hv, 'Color', 'g', 'LineWidth', 2.0, 'AutoScaleFactor', 1/4, 'ShowArrowHead','off');
                    hv = quiver( xc_plot(1:step:end), zc_plot(1:step:end),     -v21(1:step:end,1:step:end),     -v22(1:step:end,1:step:end) );%,  0.25, 'g', 'AutoScaleFactor', 1/8,  'LineWidth', 3, 'ShowArrowHead','off')
                    set(hv, 'Color', 'g', 'LineWidth', 2.0, 'AutoScaleFactor', 1/4, 'ShowArrowHead','off');
                end
            end
            
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
        
        figure(91)
        subplot(311)
        imagesc(xc_plot,zc_plot,ndx), colorbar, set(gca,'ydir','normal'), axis image
        subplot(312)
        imagesc(xc_plot,zc_plot,ndz), colorbar, set(gca,'ydir','normal'), axis image
        subplot(313)
        imagesc(xc_plot,zc_plot,sxz), colorbar, set(gca,'ydir','normal'), axis image
        
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        toc
        drawnow
    end
end

end

function VizGrid = PhaseMap( filename, VizGrid )

VizGrid.ph = hdf5read(filename,'/VizGrid/compo'); VizGrid.ph = cast(VizGrid.ph, 'double');
VizGrid.ph = (reshape(VizGrid.ph,length(VizGrid.x),length(VizGrid.z))');

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

% phases = [-1; 3; 4; 6]
phases = [-1; 0; 3; 4; 6]
for j = 1:length(phases)%max(VizGrid.ph(:))-1
    
    i = phases(j);
    Vizbuf = VizGrid.ph;
    Vizbuf(VizGrid.ph == i) = 1;
    Vizbuf(VizGrid.ph ~= i) = 69;
    [c, hb] = contour(VizGrid.x_plot, VizGrid.z_plot, Vizbuf, [1]);
    set(hb, 'Color', 'w', 'LineWidth', 0.5);
    
end
end

function  [L1,L2,s11,s21,s12,s22] = Eigs_2x2( axx, axz, azx, azz )

% Eigenvalues
Tr   = axx + azz;
Det  = axx.*azz - axz.*azx;
L1   = Tr./2 + sqrt( Tr.^2/4 - Det );
L2   = Tr./2 - sqrt( Tr.^2/4 - Det );

% Extract principal strain axis 2 -- finite strain
a22 = azz - L2;
a21 = azx;
alp = (a22./a21);    % should be correct for a normal basis vector
s11 = ones(size(a21)) .* 1./sqrt(alp.^2+1) .* L2;
s12 = alp             .* s11;
% s11 = reshape(s11,ny,nx)'; s11 = s11(1:1:size(s11,1), 1:1:size(s11,2));
% s12 = reshape(s12,ny,nx)'; s12 = s12(1:1:size(s12,1), 1:1:size(s12,2));

% Extract principal strain axis 2 -- finite strain
a22 = azz - L1;
a21 = azx;
alp = (a22./a21);    % should be correct for a normal basis vector
s21 = ones(size(a21)) .* 1./sqrt(alp.^2+1) .* L1;
s22 = alp             .* s21;
% s21 = reshape(s21,1,1)'; s21 = s21(1:1:size(s21,1), 1:1:size(s21,2));
% s22 = reshape(s22,1,1)'; s22 = s22(1:1:size(s22,1), 1:1:size(s22,2));
end

function [c11,c12,c21,c22] =  AxB_2x2( a11, a12, a21, a22, b11,b12,b21,b22 )
c11 = a11.*b11 + a12.*b21;
c12 = a11.*b12 + a12.*b22;
c21 = a21.*b11 + a22.*b21;
c22 = a21.*b12 + a22.*b22;
end

function [b11,b12,b21,b22] =  Sqrt_2x2( a11, a12, a21, a22 )
Tr  = a11 + a22;
Det = a11.*a22 - a12.*a21;
s   = sqrt(Det);
t   = sqrt(Tr + 2*s);
b11 = 1./t.*(a11 + s);
b12 = 1./t.*a12;
b21 = 1./t.*a21;
b22 = 1./t.*(a22 + s);
end

function  [aixx,aizx,aixz,aizz] = Inv_2x2( axx, axz, azx, azz )
Det  = axx.*azz - axz.*azx;
aixx = 1./Det .*azz;
aixz =-1./Det .*axz;
aizx =-1./Det .*azx;
aizz = 1./Det .*axx;
end