clear, %close all

path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE/'

istart = 1
ijump  = 1;
iend   = istart;

icount = 0;
time   = 0;
print2screen = 1;
crop   = 0;
% 
% minVx = -6e5;
% maxVx = 6e5
% minVz = -6e5;
% maxVz = 6e5

for istep=istart:ijump:iend
    
    tic
    
    figCount = 0;
    icount = icount + 1;
    filename = [path,'Output',num2str(0,'%05d'),'','.gzip.h5'];
 
    params     = hdf5read(filename,'/Model/Params');
    xg_coord   = hdf5read(filename,'/Model/xg_coord');  xg_coord  = cast(xg_coord, 'double');
    zg_coord   = hdf5read(filename,'/Model/zg_coord');  zg_coord  = cast(zg_coord, 'double');
    xc_coord   = hdf5read(filename,'/Model/xc_coord');  xc_coord  = cast(xc_coord, 'double');
    zc_coord   = hdf5read(filename,'/Model/zc_coord');  zc_coord  = cast(zc_coord, 'double');
    xvz_coord  = hdf5read(filename,'/Model/xvz_coord'); xvz_coord = cast(xvz_coord, 'double');
    zvx_coord  = hdf5read(filename,'/Model/zvx_coord'); zvx_coord = cast(zvx_coord, 'double');
    VizGrid.x  = hdf5read(filename,'/VizGrid/xviz');    VizGrid.x = cast(VizGrid.x, 'double');
    VizGrid.z  = hdf5read(filename,'/VizGrid/zviz');    VizGrid.z = cast(VizGrid.z, 'double');
    
    xc_plot = xc_coord;
    zc_plot = zc_coord;
    xg_plot = xg_coord;
    zg_plot = zg_coord;
    zvx_plot = zvx_coord;
    xvz_plot = xvz_coord;
    
    nx    = params(4);
    nz    = params(5);
    ncx   = nx-1;
    ncz   = nz-1;
    time0 = time;
    time  = params(1);
    t_vec(icount) = time; 
    
    xLabel = '';
    zLabel = '';
    TimeLabel = '';

    %--------------------------------------------------
    % plot residuals 
    %--------------------------------------------------
        
        filename = [path,'Residuals',num2str(istep,'%05d'),'','.gzip.h5'];
        VxNodes.Vx = hdf5read(filename,'/VxNodes/ru');
        VzNodes.Vz = hdf5read(filename,'/VzNodes/rv');
        VxNodes.Vx = cast(VxNodes.Vx, 'double');
        VzNodes.Vz = cast(VzNodes.Vz, 'double');
        
        if print2screen == 1
            figCount = figCount +1;
            figure(figCount), clf
        else
            figure('Visible', 'Off')
        end
        
        clf
        % Vx
        subplot(1,2,1)
        Vx = reshape(VxNodes.Vx,nx,nz+1)'; 
        Vx = Vx(2:end-1,:);
        imagesc(xg_plot/1e3, zvx_plot(2:end-1)/1e3, log10(abs(Vx)) );
        shading flat,axis xy image, colorbar;
        xlabel(xLabel), ylabel(zLabel);
%         title(['Fx at' TimeLabel, ' min = ', num2str(min(Vx(:)), '%2.4e'), ' max = ', num2str(max(Vx(:)), '%2.4e')])
        if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
        if exist('minVx', 'var') caxis([minVx maxVx]); end
        
        subplot(1,2,2)
        Vz = reshape(VzNodes.Vz,params(4)+1,params(5))';  
        Vz = Vz(:,2:end-1);
        imagesc(xvz_plot(2:end-1)/1e3,zg_plot/1e3, log10(abs(Vz)) );        
        shading flat,axis xy image, colorbar;
        xlabel(xLabel), ylabel(zLabel);
%         title(['Fz' TimeLabel, ' min = ', num2str(min(Vz(:)), '%2.4e'), ' max = ', num2str(max(Vz(:)), '%2.4e')])
         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
        if exist('minVz', 'var') caxis([minVz maxVz]); end
        
        
%         eta_s = hdf5read(filename,'/Vertices/eta'); eta_s = cast(eta_s, 'double');
%         eta_n = hdf5read(filename,'/Centers/eta');  eta_n = cast(eta_n, 'double');
% 
%         eta_s = reshape(eta_s,params(4),params(5))';
%         eta_n = reshape(eta_n,params(4)-1,params(5)-1)';
%         
%         rho_s = hdf5read(filename,'/Vertices/rho'); rho_s = cast(rho_s, 'double');
%         rho_n = hdf5read(filename,'/Centers/rho');  rho_n = cast(rho_n, 'double');
%         rho_s = reshape(rho_s,params(4),params(5))';
%         rho_n = reshape(rho_n,params(4)-1,params(5)-1)';
% 
%         % Nodes
%         if print2screen == 1
%             figCount = figCount +1;
%             f=figure(figCount), clf
%         else
%             figure('Visible', 'Off')
%         end
%          
%         subplot(211)
%         hold on
%         imagesc(xg_plot, zg_plot, log10(eta_s));
%         shading flat, axis xy image, colorbar;
%         title(['Viscosity (log10) at' TimeLabel])
%         xlabel(xLabel), ylabel(zLabel);
%         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%         if exist('min_eta', 'var') caxis([min_eta max_eta]); end
%                  
%         subplot(212)
%         hold on
%         imagesc(xc_plot, zc_plot, log10(eta_n));
%         shading flat, axis xy image, colorbar;
%         title(['Viscosity (log10) at' TimeLabel])
% 
%         xlabel(xLabel), ylabel(zLabel);
%         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
        
%         % Nodes
%         if print2screen == 1
%             figCount = figCount +1;
%             f=figure(figCount), clf
%         else
%             figure('Visible', 'Off')
%         end
%          
%         subplot(211)
%         hold on
%         imagesc(xg_plot, zg_plot, (rho_s));
%         shading flat, axis xy image, colorbar;
%         title(['Density s at' TimeLabel])
%         xlabel(xLabel), ylabel(zLabel);
%         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
%         if exist('min_eta', 'var') caxis([min_eta max_eta]); end
%                  
%         subplot(212)
%         hold on
%         imagesc(xc_plot, zc_plot, (rho_n));
%         shading flat, axis xy image, colorbar;
%         title(['Density n at' TimeLabel])
% 
%         xlabel(xLabel), ylabel(zLabel);
%         if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end

norm(Vx)
norm(Vz)   
    end