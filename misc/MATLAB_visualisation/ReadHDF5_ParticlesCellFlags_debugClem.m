clear
% close all
clc

path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE/'
mdoodz6       = 1;

cd(path)

% File numbers
istart = 0;
ijump  = 1;
iend   = 0;

iter   = 0;

% Visualisation options
MarkSize      = 5e0;
printfig      = 0;
print2screen  = 1;
 
format        = '-dpng';
res           = '-r200';
color_map     = 'piston0';
file_suffix   = '';

particles = 1;
topo      = 0;
vel_plot  = 0;
tags      = 0;
residuals = 0;
numbering = 0;

% Colorbar limits
minP   = -2e5;


maxP   = 1e5;
minRho = 2400;
maxRho = 3100;
maxS   = 3;
minS   = -4;

minTopo = -5e3;
maxTopo = 3e3;

% Size of the window
crop   = 0;

xmin   =  -4e3;
xmax   =  4e3;
zmin   = -4e3;
zmax   =  1e3;

% Arrays
gamma  = zeros(10000,1);
error  = zeros(10000,1);
etapp  = zeros(10000,1);
icount = 0;

figCount = 0;

% Create directories
if (exist('./Fig_Phases', 'dir') == 0)
    mkdir ./Fig_Phases
end


%% LOOP ON OUTPUT FILES

for istep=istart:ijump:iend
    
    icount = icount + 1;  
    
    filename = [path,'Output',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];
    filename_p = [path,'Particles',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];
    
%      filename = [path,'Output_BeforeSolve',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];
%     filename_p = [path,'Particles_BeforeSolve',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];


%     filename = [path,'Outputxx',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];
%     filename_p = [path,'Particlesxx',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];

%     
%     filename = [path,'DEBUG',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];
%     filename_p = [path,'DEBUGParticles',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];

%     filename = [path,'Output',num2str(istep-1,'%05d'),file_suffix,'.gzip.h5'];
%     filename_p = [path,'AfterReseeding',num2str(istep,'%05d'),file_suffix,'.gzip.h5'];


    fileinfo = hdf5info(filename);
    
    params     = hdf5read(filename,'/Model/Params');
    
    xg_coord   = hdf5read(filename,'/Model/xg_coord');  xg_coord  = cast(xg_coord, 'double');
    zg_coord   = hdf5read(filename,'/Model/zg_coord');  zg_coord  = cast(zg_coord, 'double');
    xc_coord   = hdf5read(filename,'/Model/xc_coord');  xc_coord  = cast(xc_coord, 'double');
    zc_coord   = hdf5read(filename,'/Model/zc_coord');  zc_coord  = cast(zc_coord, 'double');
    xvz_coord  = hdf5read(filename,'/Model/xvz_coord'); xvz_coord = cast(xvz_coord, 'double');
    zvx_coord  = hdf5read(filename,'/Model/zvx_coord'); zvx_coord = cast(zvx_coord, 'double');
    
    
    xc_plot = xc_coord;
    zc_plot = zc_coord;
    xg_plot = xg_coord;
    zg_plot = zg_coord;
    zvx_plot = zvx_coord;
    xvz_plot = xvz_coord;
        
    nx = params(4);
    nz = params(5);
    ncx = nx-1;
    ncz = nz-1;
    dx = xg_coord(2) - xg_coord(1);
    dz = zg_coord(2) - zg_coord(1);
    time = params(1);
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
    
    
    params(1) = params(1) * 365 * 24 * 3600;
    
    %--------------------------------------------------
    % plot particules + phases
    %--------------------------------------------------
    
    if particles
       
        Part.x     = hdf5read(filename_p,'/Particles/x');
        Part.z     = hdf5read(filename_p,'/Particles/z');
        Part.ph    = hdf5read(filename_p,'/Particles/phase');
        Part.gen    = hdf5read(filename_p,'/Particles/generation');
        
%         zmin_RTFS(icount) = min(Part.z(Part.ph==0));
%         t_vec(icount) = params(1) / (365 * 24 * 3600);
        
        if topo
            xtopo = hdf5read(filename_p,'/Topo/x');
            ztopo = hdf5read(filename_p,'/Topo/z');
            phase = hdf5read(filename_p,'/Topo/phase');
            vxtopo = hdf5read(filename_p,'/Topo/vx');
            vztopo = hdf5read(filename_p,'/Topo/vz');
            xtopo = cast(xtopo, 'double');
            ztopo = cast(ztopo, 'double');
            height = hdf5read(filename_p,'/Topo/height');
            vxs = hdf5read(filename_p,'/Topo/vxsurf');
            vzs = hdf5read(filename_p,'/Topo/vzsurf');
            
            xtopo_ini = hdf5read(filename_p,'/Topo_ini/x');
            ztopo_ini = hdf5read(filename_p,'/Topo_ini/z');
            phase_ini = hdf5read(filename_p,'/Topo_ini/phase');
            vxtopo_ini = hdf5read(filename_p,'/Topo_ini/vx');
            vztopo_ini = hdf5read(filename_p,'/Topo_ini/vz');
            xtopo_ini = cast(xtopo_ini, 'double');
            ztopo_ini = cast(ztopo_ini, 'double');
            height_ini = hdf5read(filename_p,'/Topo_ini/height');
            vxs_ini = hdf5read(filename_p,'/Topo_ini/vxsurf');
            vzs_ini = hdf5read(filename_p,'/Topo_ini/vzsurf');

        end

        if print2screen == 1
            figure(1), clf
        else
            figure('Visible', 'Off')
        end
        
%         if topo ==1
%         subplot (212)
%         end
        
        hold on
                mesh(xg_coord./1e3, zg_coord./1e3,zeros(nz,nx))
%                 mesh(xc_coord./1e3, zc_coord./1e3,zeros(ncz,ncx))
% %         [x2,z2]=meshgrid(xc_coord(:),zc_coord(:));
%         plot(x2./1e3, z2./1e3,'x')

        
        Part.x = Part.x/1e3;
        Part.z = Part.z/1e3;
        
        if crop == 1
            ind = Part.x>xmin & Part.x<xmax & Part.z<zmax & Part.z>zmin;
            Part.x  = Part.x(ind);
            Part.z  = Part.z(ind);
            Part.ph = Part.ph(ind);
            Part.gen= Part.gen(ind);
        end
        
        if strcmp(color_map, 'inc2boules') == 1
            hold on, plot(Part.x(Part.ph==0),Part.z(Part.ph==0),'.','MarkerSize',MarkSize,'color',[0.028 0.54 0.295]);
            hold on, plot(Part.x(Part.ph==3),Part.z(Part.ph==3),'k.','MarkerSize',MarkSize);
            hold on, plot(Part.x(Part.ph==1),Part.z(Part.ph==1),'r.','MarkerSize',MarkSize);
            hold on, plot(Part.x(Part.ph==2),Part.z(Part.ph==2),'w.','MarkerSize',MarkSize);
        elseif strcmp(color_map, 'piston0') == 1
            hold on, plot(Part.x(Part.ph==0),Part.z(Part.ph==0),'.','MarkerSize',MarkSize,'color',[56/255 117/255 215/255]);
            hold on, plot(Part.x(Part.ph==1),Part.z(Part.ph==1),'.','MarkerSize',MarkSize,'color',[1/255 31/255 113 /255]);
            hold on, plot(Part.x(Part.ph==2),Part.z(Part.ph==2),'r.','MarkerSize',MarkSize);
            hold on, plot(Part.x(Part.ph==3),Part.z(Part.ph==3),'.','MarkerSize',MarkSize,'color',[248/255 216/255 71/255] );
            hold on, plot(Part.x(Part.ph==4),Part.z(Part.ph==4),'k.','MarkerSize',MarkSize,'color',[237/255 237/255 237/255]);
            hold on, plot(Part.x(Part.ph==5),Part.z(Part.ph==5),'w.','MarkerSize',MarkSize); % Fp
            hold on, plot(Part.x(Part.ph==6),Part.z(Part.ph==6),'.','MarkerSize',MarkSize,'color',[237/255 237/255 237/255] ); % Pl
            hold on, plot(Part.x(Part.ph==7),Part.z(Part.ph==7),'k.','MarkerSize',MarkSize); % Bt
            hold on, plot(Part.x(Part.ph==8),Part.z(Part.ph==8),'g.','MarkerSize',MarkSize); % Bt
        elseif strcmp(color_map, 'piston') == 1
                 map = [115/255     115/255     115/255;...
                200/255     200/255     200/255;...
                130/255     0/255     10/255;...
                245/255     0/255     0/255;...
                75/255     75/255     75/255;....
                155/255     155/255     155/255;...
                255/255     255/255     255/255;...
                1/255     0/255     0/255;];
            for i=1:size(map,1)
                 hold on,
                 ind0  = Part.ph==i-1 & Part.gen==0;
                 ind1  = Part.ph==i-1 & Part.gen==1;
                 sum(ind1)
                 sum(Part.gen)
                 plot(Part.x(ind0)./1e0,Part.z(ind0)./1e0,'.','MarkerSize',MarkSize,'color',map(i,:));
                 plot(Part.x(Part.gen==1)./1e3,Part.z(Part.gen==1)./1e0,'*','MarkerSize',2*MarkSize,'color',map(i,:));
            end
            plot(Part.x(Part.ph==-1)./1e0,Part.z(Part.ph==-1)./1e0,'g.','MarkerSize',MarkSize);
        end
        if topo
            ind = phase==0;
            plot( xtopo(ind)./1e3, ztopo(ind)./1e3, 'dr' )
            ind = phase==-1;
            plot( xtopo(ind)./1e3, ztopo(ind)./1e3, 'xr' )
            % advected initial topo
            ind = phase_ini==0;
            plot( xtopo_ini(ind)./1e3, ztopo_ini(ind)./1e3, 'dk' )
            ind = phase_ini==-1;
            plot( xtopo_ini(ind)./1e3, ztopo_ini(ind)./1e3, 'xk' )
        end
       
%         set(gca,'position',[0 0 1 1],'units','normalized')
%         axis off
           axis image
           title(['Phases at' TimeLabel]);
             if crop == 1
            xlim([xmin xmax])
            ylim([zmin zmax])
        end
           drawnow
%            if topo
%            subplot(211)
%            
%             plot( xtopo./1e3, (ztopo-mean( ztopo(:)))./1e0, '.-', xg_coord./1e3, (height-mean(height(:)))./1e0, 'or'  )
%             title('Topography [km]')
% %             axis([xmin/1e3 xmax/1e3 minTopo/1e3 maxTopo/1e3 ])
%            end

        
       
        
%         if printfig == 1
%             if crop == 1
%                 print([path, './Fig_Phases/Crop_PhasesPart', num2str(istep,'%05d'),file_suffix], '-dtiff', res)
%             else
%                 print([path, './Fig_Phases/Fig_PhasesPart', num2str(istep,'%05d'),file_suffix], '-dtiff', res)
%             end
%             close all
%         end


        if topo
            figure(2)
            subplot(311), hold on
            plot( xtopo, vxtopo, '.r', xg_coord, vxs, 'or' )
            plot( xtopo_ini, vxtopo_ini, '.k', xg_coord, vxs_ini, 'ok' )
            title('Vx')

            subplot(312), hold on
            plot( xtopo, vztopo, '.r', xvz_coord, vzs, 'or'  )
            plot( xtopo_ini, vztopo_ini, '.k', xvz_coord, vzs_ini, 'ok'  )
            title('Vz')

            subplot(313), hold on
            plot( xtopo, ztopo, '.r', xg_coord, height, 'or'  )
            plot( xtopo_ini, ztopo_ini, 'sg', xg_coord, height_ini, 'ok'  )
            title(['mean(h) = ', num2str(mean(height))])
        end


    end
        
        
        
%         if tags
%             Tagp = hdf5read(filename,'/Centers/BCp'); % Tagp = cast(Tagp, 'double');
%             Tagp = reshape(Tagp,nx-1,nz-1)';
%             Tagp1 = Tagp(:);
%             [x2,z2] = meshgrid(xc_coord,zc_coord);
%             x1    = x2(:);
%             z1    = z2(:);
%             figure(3), clf
%             ind0 = Tagp1 == -1;
% %             imagesc(xc_coord/1e3,zc_coord/1e3,Tagp)
% %                         plot(x1./1e3,z1./1e3,'.','MarkerSize',MarkSize,'color','g');
%             plot(x1(ind0)./1e3,z1(ind0)./1e3,'.','MarkerSize',MarkSize,'color','k');
%             axis image
%             colorbar
% % %             caxis([-1 2])
%         end
% %         

        if tags == 1
            
            tag_s = hdf5read(filename,'/Flags/tag_s');
            tag_n = hdf5read(filename,'/Flags/tag_n');
            tag_u = hdf5read(filename,'/Flags/tag_u');
            tag_v = hdf5read(filename,'/Flags/tag_v');

            if print2screen == 1
                figure(1), clf
            else
                figure('Visible', 'Off')
            end
           
%             subplot(211)
%             xtopoc = hdf5read(filename,'/Topo/x_mark');
%             Vxtopo = hdf5read(filename,'/Topo/Vx_mark');
%             Vztopo = hdf5read(filename,'/Topo/Vz_mark');
%             plot( xtopo, Vxtopo, '*m' )
%                             xlim([xmin xmax])
%                 ylim([-5e-13 5e-13])
%             
%             subplot(212), 
            hold on
            mesh(xg_coord, zg_coord, zeros(length(zg_coord),length(xg_coord)), 'edgecolor', [0.7 0.7 0.7] )
            % axis image

            [xc2, zc2] = ndgrid(xc_coord, zc_coord);
            xc2 = xc2(:);
            zc2 = zc2(:);
            tag_n = tag_n(:);
            plot(xc2(tag_n==-1), zc2(tag_n==-1), 'or')
            plot(xc2(tag_n==31), zc2(tag_n==31), 'og')


            [xn2, zvx2] = ndgrid(xg_coord, zvx_coord);
            xn2 = xn2(:);
            zvx2 = zvx2(:);
            tag_u = tag_u(:);
            plot(xn2(tag_u==-1), zvx2(tag_u==-1), '+b', 'MarkerSize', 5)
            plot(xn2(tag_u==11), zvx2(tag_u==11), '+b', 'MarkerSize', 5)
            plot(xn2(tag_u==0), zvx2(tag_u==0), '+b', 'MarkerSize', 5)

            [xvz2, zn2] = ndgrid(xvz_coord, zg_coord);
            xvz2 = xvz2(:);
            zn2 = zn2(:);
            tag_v = tag_v(:);
            plot(xvz2(tag_v==-1), zn2(tag_v==-1), 'dk')
            plot(xvz2(tag_v==11), zn2(tag_v==11), 'dk')
            plot(xvz2(tag_v==0), zn2(tag_v==0), 'dk')

            [xn2, zn2] = ndgrid(xg_coord, zg_coord);

            plot(xn2(tag_s==-1), zn2(tag_s==-1), '.m')
            % plot(xn2(tag_s==30), zn2(tag_s==30), '.y')
            plot(xc2(tag_n==0), zc2(tag_n==0), 'dr')

            if topo == 1
                if mdoodz6 == 1
                    xtopo = hdf5read(filename,'/Topo/x_mark');
                    ztopo = hdf5read(filename,'/Topo/z_mark');
                else
                    xtopo = hdf5read(filename,'/Topo/x');
                    ztopo = hdf5read(filename,'/Topo/z');
                end
                plot( xtopo, ztopo, '*m' )
            end
            if crop == 1
                xlim([xmin xmax])
                ylim([zmin zmax])
            end
            
            
            if printfig == 1
                print([path, './Fig_Phases/Fig_DEBUG', num2str(istep,'%05d'),file_suffix], '-dtiff', res)
                close all
            end

            
        end
        
        if residuals == 1
            
            %--------------------------------------------------
            % plot residuals
            %--------------------------------------------------
            
            filename = [path,'Residuals',num2str(iter,'%05d'),'','.gzip.h5'];
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
            xlabel('x'), ylabel('z');
            %         title(['Fx at' TimeLabel, ' min = ', num2str(min(Vx(:)), '%2.4e'), ' max = ', num2str(max(Vx(:)), '%2.4e')])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minVx', 'var') caxis([minVx maxVx]); end
            
            subplot(1,2,2)
            Vz = reshape(VzNodes.Vz,params(4)+1,params(5))';
            Vz = Vz(:,2:end-1);
            imagesc(xvz_plot(2:end-1)/1e3,zg_plot/1e3, log10(abs(Vz)) );
            shading flat,axis xy image, colorbar;
             xlabel('x'), ylabel('z');
            %         title(['Fz' TimeLabel, ' min = ', num2str(min(Vz(:)), '%2.4e'), ' max = ', num2str(max(Vz(:)), '%2.4e')])
            if crop == 1 xlim([lim.xmin lim.xmax]); ylim([lim.zmin lim.zmax]); end
            if exist('minVz', 'var') caxis([minVz maxVz]); end
            
       
            
            norm(Vx)
            norm(Vz)
        end
        
        if numbering == 1
            filename = [path,'MatrixDecoupled.gzip_1cpu.h5'];
            eqnu = hdf5read(filename,'/numbering/eqn_u');
            eqnv = hdf5read(filename,'/numbering/eqn_v');
            eqnp = hdf5read(filename,'/numbering/eqn_p');
            eqnu = reshape(eqnu,nx,nz+1)';
            eqnv = reshape(eqnv,nx+1,nz)';
            eqnp = reshape(eqnp,nx-1,nz-1)';
            
            if print2screen == 1
                figCount = figCount +1;
                figure(figCount), clf
            else
                figure('Visible', 'Off')
            end
            
            subplot(211)
            imagesc(xg_plot/1e3, zvx_plot(2:end-1)/1e3, eqnu(2:end-1,:))
            shading flat,axis xy image, colorbar;
            xlabel('x'), ylabel('z');
            
            subplot(212)
            imagesc(xc_plot/1e3, zc_plot/1e3, eqnp)
            shading flat,axis xy image, colorbar;
            xlabel('x'), ylabel('z');
            
        end
            
        
        clear Part.x Part.z Part.ph Part
        colormap(jet)
        drawnow
        
end


