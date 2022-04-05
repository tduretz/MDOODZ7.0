clear 
close all
colormap(jet);
clc

path = '/Users/tduretz/REPO_GIT/MDOODZ6.0/SOURCE/';

cd(path)

% File
istart = 1;
ijump  = 1;
iend   = 25;
icount = 0;

%--------------------------------------------------
% what do dou you want to plot:
%--------------------------------------------------

pre_plot        = 1;

% Visualisation options
Ftsz          = 15;
file_suffix   = '';

%% LOOP ON OUTPUT FILES
for istep=istart:ijump:iend
    
    tic
    
    figCount = 0;
    icount = icount + 1;
    
    filename = [path,'Output',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];
    
    params     = hdf5read(filename,'/Model/Params');
    xg_coord   = hdf5read(filename,'/Model/xg_coord');  xg_coord  = cast(xg_coord, 'double');
    zg_coord   = hdf5read(filename,'/Model/zg_coord');  zg_coord  = cast(zg_coord, 'double');
    xc_coord   = hdf5read(filename,'/Model/xc_coord');  xc_coord  = cast(xc_coord, 'double');
    zc_coord   = hdf5read(filename,'/Model/zc_coord');  zc_coord  = cast(zc_coord, 'double');
    xvz_coord  = hdf5read(filename,'/Model/xvz_coord'); xvz_coord = cast(xvz_coord, 'double');
    zvx_coord  = hdf5read(filename,'/Model/zvx_coord'); zvx_coord = cast(zvx_coord, 'double');
    VizGrid.x  = hdf5read(filename,'/VizGrid/xviz');    VizGrid.x = cast(VizGrid.x, 'double');
    VizGrid.z  = hdf5read(filename,'/VizGrid/zviz');    VizGrid.z = cast(VizGrid.z, 'double');

    nx    = params(4);
    nz    = params(5);
    ncx   = nx-1;
    ncz   = nz-1;
    time  = params(1);
    dx    = xg_coord(2) - xg_coord(1);
    dz    = zg_coord(2) - zg_coord(1);

    %--------------------------------------------------
    % plot pressure
    %--------------------------------------------------
    if (pre_plot==1)
        
        % Extract fields
        P    = hdf5read(filename,'/Centers/P');    P    = cast(   P , 'double');
        P    = reshape(   P,params(4)-1,params(5)-1)';
        sxxd = hdf5read(filename,'/Centers/sxxd'); sxxd = cast(sxxd , 'double');
        sxxd = reshape(sxxd,params(4)-1,params(5)-1)';
        Vx   = hdf5read(filename,'/VxNodes/Vx'); Vx = cast(Vx, 'double');
        Vx   = reshape(Vx,nx,nz+1)';
        etac = hdf5read(filename,'/Centers/eta_n');  etac = cast(etac, 'double');
        etac = reshape(etac,nx-1,nz-1)';
       
        % Compute total strain rate (including divergence)
        exx  = diff(Vx,1,2)/dx;
        
        % Plot
        figure(10),
        eta_layer = 1.0;
        beta      = 1e5;
        tmax      = eta_layer*beta;
        dmus      = etac(fix(nz/2),1) - etac(1,1);
        ffs       = abs( 2*exx(fix(nz/2),1)* dmus);
        dP        = P(fix(nz/2),1) - P(1,1);   
        dTxx      = sxxd(fix(nz/2),1) - sxxd(1,1);
        dPa       =  ffs*(1 - exp(-3/4*time/tmax));
        dTxxa     = -ffs*( 1 - 1/2*exp(-3/4*time/tmax));
        subplot(211), hold on
        plot( time/tmax,  dTxx /ffs, 'ok' )
        plot( time/tmax,  dTxxa/ffs, 'xr' )
        xlabel('$t$','interpreter','latex'); ylabel('$d\tau_{xx}$','interpreter','latex'); l1=legend('$\tau{xx}$','$\tau{xx}$ anal.'); set(l1,'interpreter','latex', 'box', 'off')
        axis([0 5 -1 -0.5])
        set(gca, 'FontSize', 15)
        subplot(212),  hold on
        plot( time/tmax,  dP/ffs, 'ok' )
        plot( time/tmax, dPa/ffs, 'xr' )
        xlabel('$t$','interpreter','latex'); ylabel('$dP$','interpreter','latex'); l2=legend('$dP$', '$dP$ anal.'); set(l2,'interpreter','latex','location','southeast', 'box', 'off')
        axis([0 5 0 1])
        set(gca, 'FontSize', Ftsz)

    end
    
    drawnow
    toc
end