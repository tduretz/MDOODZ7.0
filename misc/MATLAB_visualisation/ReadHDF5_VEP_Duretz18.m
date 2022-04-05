
function Main

clear
clear all
% close all
clc

DEBUG = 0;

MarkSize=1e0 ;


path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE//'
M2Di_EP = load('./BENCHMARKS_data/DataM2Di_EP_test01');


cd(path)

% File
istart = 100;
ijump  = 1
iend   = 100;

%--------------------------------------------------
% what do you want to plot:
%--------------------------------------------------
vel_divergence  = 0;
pre_plot        = 0;
stress_inv      = 1;
stress_evol     = 0;
divergence      = 0;

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
MaskAir       = 1;
ColorFabio    = 1;
y             = 365.25*24*3600;
My            = 1e6*y;
topo_ref      = 0; %6370
mdoodz6       = 1;

% Contour
eps_val = 5e-15;

minRho = 3000;
maxRho = 3600;

% Colorbar limits
minP = -6e6;
maxP =  4e6;

minEii = -16;
maxEii = -13;

% Size of the window
crop       = 0;

lim.xmin   = -175;
lim.xmax   =  175;
lim.zmin   = -175;
lim.zmax   =  175;

% Arrays
t_vec  = zeros(10000,1);
icount = 0;
time   = 0;

short       = zeros( 10000,1);
siivec      = zeros( 10000,1);
eiivec      = zeros( 10000,1);
eiimaxvec   = zeros( 10000,1);

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
        % plot velocity divergence
        %--------------------------------------------------
        if ( vel_divergence == 1 )
            
            VxNodes.Vx = hdf5read(filename,'/VxNodes/Vx'); VxNodes.Vx = cast(VxNodes.Vx, 'double');
            VzNodes.Vz = hdf5read(filename,'/VzNodes/Vz'); VzNodes.Vz = cast(VzNodes.Vz, 'double');
            Vx = reshape(VxNodes.Vx,nx,nz+1)';
            Vz = reshape(VzNodes.Vz,nx+1,nz)';
            
            % Velocity divergence
            dVxdx = (Vx(:,2:end) - Vx(:,1:end-1)) ./ (xc_coord(2)-xc_coord(1));
            dVzdz = (Vz(2:end,:) - Vz(1:end-1,:)) ./ (zc_coord(2)-zc_coord(1));
            divV  = abs(dVxdx(2:end-1,:) + dVzdz(:,2:end-1));
            
            if crop == 1
                [divV, xc_plot, zc_plot] = CropCellArray( divV, xc_plot, zc_plot, lim );
            end
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            imagesc( xc_plot, zc_plot, log10(divV) )
            shading flat,axis xy image, colorbar;
            xlabel('X'), ylabel('Z')
            title(['Velocity divergence at' TimeLabel])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            
            if printfig == 1
                print([path, './Fig_VelocityDivergence/Fig_VelocityDivergence', num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
        end
        
        
        %--------------------------------------------------
        % plot pressure
        %--------------------------------------------------
        if (pre_plot==1)
            Cent.P  = hdf5read(filename,'/Centers/P');
            Cent.P  = cast(Cent.P , 'double');
            P = reshape(Cent.P,params(4)-1,params(5)-1)';

            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            hold on
            imagesc( xc_plot, zc_plot, P)
            shading flat,axis xy image, colorbar;
            hold off

            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minP', 'var') caxis([minP maxP]); end
            xlabel(xLabel), ylabel(zLabel);
            title(['min P =', num2str(min(P(:))), ' max P = ', num2str(max(P(:)))])
            
            if printfig == 1
                print([path,'./Fig_Pressure/Fig_Pressure',num2str(istep,'%05d'),file_suffix], format, res)
                close all
            end
            
        end

        %--------------------------------------------------
        % plot stress and strain rate invariants evolution
        %--------------------------------------------------
        if ( stress_evol == 1 )
            
            Cent.sxxd  = hdf5read(filename,'/Centers/sxxd'); Cent.sxxd = cast(Cent.sxxd, 'double');
            Cent.P     = hdf5read(filename,'/Centers/P');    Cent.P    = cast(Cent.P, 'double');
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
            
            sxxd = (reshape(Cent.sxxd,params(4)-1,params(5)-1)');
            P    = (reshape(Cent.P,params(4)-1,params(5)-1)');
            szzd = -sxxd;
            sxz  = (reshape(Vert.sxz, params(4)  ,params(5)  )');
            sII  = sqrt( 0.5*(2*sxxd.^2 + 0.5*(sxz(1:end-1,1:end-1).^2 + sxz(2:end,1:end-1).^2 + sxz(1:end-1,2:end).^2 + sxz(2:end,2:end).^2 ) ) );
            
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            HSc = 2*sxzc.*exzc + sxxd.*exx + szzd.*ezz;
            
            sII = sII(sII>1e-8);
            vol = length(sII)*dx*dz;
            stress = sum(sII(:).*(dx*dz))/vol;
            if istep == 0
                stress = 0;
            end
            
            mean_stress(icount) = mean(sII(:));
            min_stress(icount)  = min(sII(:));
            max_stress(icount)  = max(sII(:));
            
            mean_press(icount) = mean(P(:));
            min_press(icount)  = min(P(:));
            max_press(icount)  = max(P(:));

            figure(20),clf,colormap('jet'),set(gcf,'Color','white')
            subplot(211), hold on
            plot((1:icount)*dt, (mean_stress(1:icount)),'b-')
            plot((1:icount)*dt, (max_stress(1:icount)),'r-')
            plot((1:icount)*dt, (min_stress(1:icount)),'g-')
            plot((1:icount)*dt, M2Di_EP.Tiivec(1:icount),'ob')
            plot((1:icount)*dt, M2Di_EP.Tiimin(1:icount),'og')
            plot((1:icount)*dt, M2Di_EP.Tiimax(1:icount),'or')
            xlabel('Time'),ylabel('Stress')
            subplot(212), hold on
            plot((1:icount)*dt, (mean_press(1:icount)),'b-')
            plot((1:icount)*dt, (max_press(1:icount)),'r-')
            plot((1:icount)*dt, (min_press(1:icount)),'g-')
            plot((1:icount)*dt, M2Di_EP.Pvec(1:icount),'bo')
            plot((1:icount)*dt, M2Di_EP.Pminv(1:icount),'go')
            plot((1:icount)*dt, M2Di_EP.Pmaxv(1:icount),'ro')
            xlabel('Time'),ylabel('Pressure')
            drawnow
            
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
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %--------------------------------------------------
        % plot stress and strain rate invariants
        %--------------------------------------------------
        if ( stress_inv == 1 )
            
            sxxd = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd, 'double');
            szzd = hdf5read(filename,'/Centers/szzd'); szzd = cast(szzd, 'double');
            sxz  = hdf5read(filename,'/Vertices/sxz'); sxz  = cast(sxz,  'double');
            exxd = hdf5read(filename,'/Centers/exxd'); exxd = cast(exxd, 'double');
            ezzd = hdf5read(filename,'/Centers/ezzd'); ezzd = cast(ezzd, 'double');
            exz  = hdf5read(filename,'/Vertices/exz'); exz  = cast(exz,  'double');
            
            exxd = reshape(exxd,params(4)-1,params(5)-1)';
            ezzd = reshape(ezzd,params(4)-1,params(5)-1)';
            eyyd = -(exxd + ezzd);
            exz  = (reshape(exz, params(4)  ,params(5)  )');
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            eII  = sqrt( 1/2*(exxd.^2 + eyyd.^2 + ezzd.^2) + exzc.^2);
            
            sxxd = (reshape(sxxd,params(4)-1,params(5)-1)');
            szzd = (reshape(szzd,params(4)-1,params(5)-1)');
            syyd = -(sxxd + szzd);
            sxz  = (reshape(sxz, params(4)  ,params(5)  )');
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            sII  = sqrt( 1/2*(sxxd.^2 + syyd.^2 + szzd.^2) + sxzc.^2);
            
            exzc = 0.25*(exz(1:end-1,1:end-1) + exz(2:end,1:end-1) + exz(1:end-1,2:end) + exz(2:end,2:end));
            sxzc = 0.25*(sxz(1:end-1,1:end-1) + sxz(2:end,1:end-1) + sxz(1:end-1,2:end) + sxz(2:end,2:end));
            Hs   = 2*sxzc.*exzc + sxxd.*exxd + szzd.*ezzd + syyd.*eyyd;
            
            load('roma.mat')
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
                        
            % eII
            if ( H>L )
                subplot(1,3,1)
            else
                subplot(3,1,1)
            end
%             colormap(flipud(roma));
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
            if exist('minEii', 'var') caxis([(minEii) (maxEii)]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
            
            % sII
            if ( H>L )
                subplot(1,3,2)
            else
                subplot(3,1,2)
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
            
            % Hs
            if ( H>L )
                subplot(1,3,3)
            else
                subplot(3,1,3)
            end
            imagesc( xc_plot, zc_plot, Hs )
            shading flat,axis xy image, colorbar;
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, Hs, 'visible', 'off' )
            hold off
            xlabel(xLabel), ylabel(zLabel);
            title(['Hs at' TimeLabel, ' min = ', num2str((min(Hs(:))), '%2.2e' ), ' Pa/s', ' max = ', num2str((max(Hs(:))), '%2.2e' ), ' Pa/s' ])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%             if exist('minSii', 'var') caxis([(minSii) (maxSii)]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/Fig_Stress_StrainRate',num2str(istep,'%05d'),file_suffix], format, res)
                close all
                
            end
        end
        
        %--------------------------------------------------
        % plot stress and strain rate invariants
        %--------------------------------------------------
        if ( divergence == 1 )
            
            div_tot = hdf5read(filename,'/Centers/divu'   ); div_tot = cast(div_tot, 'double'); div_tot = reshape(div_tot,params(4)-1,params(5)-1)';
            div_el  = hdf5read(filename,'/Centers/divu_el'); div_el  = cast(div_el , 'double'); div_el  = reshape(div_el ,params(4)-1,params(5)-1)';
            div_pl  = hdf5read(filename,'/Centers/divu_pl'); div_pl  = cast(div_pl , 'double'); div_pl  = reshape(div_pl ,params(4)-1,params(5)-1)';
            div_th  = hdf5read(filename,'/Centers/divu_th'); div_th  = cast(div_th , 'double'); div_th  = reshape(div_th ,params(4)-1,params(5)-1)';

            load('roma.mat')
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
                        
            % div_tot
            if ( H>L )
                subplot(2,2,1)
            else
                subplot(2,2,1)
            end
%             colormap(flipud(roma));
            imagesc( xc_plot, zc_plot, div_tot )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, div_tot, 'visible', 'off' )
            hold off
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['div_{tot} at' TimeLabel, ' min = ', num2str((min(div_tot(:))), '%2.2e'  ), ' s^{-1}', ' max = ', num2str((max(div_tot(:))), '%2.2e'  ), ' s^{-1}' ])
            caxis([-0.25e-14 0.25e-14])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
            
            % div_el
            if ( H>L )
                subplot(2,2,2)
            else
                subplot(2,2,2)
            end
%             colormap(flipud(roma));
            imagesc( xc_plot, zc_plot, div_el  )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, div_el , 'visible', 'off' )
            hold off
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['div_{el}  at' TimeLabel, ' min = ', num2str((min(div_el(:))), '%2.2e'  ), ' s^{-1}', ' max = ', num2str((max(div_el(:))), '%2.2e'  ), ' s^{-1}' ])
            caxis([-0.25e-14 0.25e-14])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
            
            % div_pl
            if ( H>L )
                subplot(2,2,3)
            else
                subplot(2,2,3)
            end
%             colormap(flipud(roma));
            imagesc( xc_plot, zc_plot, div_pl  )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, div_pl , 'visible', 'off' )
            hold off
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['div_{pl}  at' TimeLabel, ' min = ', num2str((min(div_pl(:))), '%2.2e'  ), ' s^{-1}', ' max = ', num2str((max(div_pl(:))), '%2.2e'  ), ' s^{-1}' ])
            caxis([-0.25e-14 0.25e-14])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
            
            % div_th
            if ( H>L )
                subplot(2,2,4)
            else
                subplot(2,2,4)
            end
%             colormap(flipud(roma));
            imagesc( xc_plot, zc_plot, div_th  )
            hold on
            if Ccontours == 1; AddCompoContours( filename, VizGrid, crop, lim  ); end
            imagesc( xc_plot, zc_plot, div_th , 'visible', 'off' )
            hold off
            shading flat,axis xy image, colorbar;
            xlabel(xLabel), ylabel(zLabel);
            title(['div_{th}  at' TimeLabel, ' min = ', num2str((min(div_th(:))), '%2.2e'  ), ' s^{-1}', ' max = ', num2str((max(div_th(:))), '%2.2e'  ), ' s^{-1}' ])
            caxis([-0.25e-14 0.25e-14])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if MaskAir==1, patch(x_tab(:,id)', z_tab(:,id)', repmat(f,1,4)', 'EdgeColor', 'none','FaceColor','w' ); end
            axis([min(xg_plot) max(xg_plot) min(zg_plot) max(zg_plot)])
            
            figure(22),
            imagesc( xc_plot, zc_plot, div_tot - div_el - div_pl - div_th ), shading flat,axis xy image, colorbar;
            
            if printfig == 1
                print([path,'./Fig_Stress_StrainRate/Fig_Divergences',num2str(istep,'%05d'),file_suffix], format, res)
                close all
                
            end
        end
        
        toc
        drawnow
    end
end

% print('PressureNoMarkers', '-dpng', '-r300')

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