clear
close all

%**** PATH ****%
path = '/Users/tduretz/REPO_GIT/MDOODZ6.0/SOURCE/'

phase = 0; 

%**** READ ****%
filename   = [path,'DefMap', num2str(phase,'%02d'),'.gzip.h5'];
params     = hdf5read(filename,'/model/params');
T          = hdf5read(filename,'/arrays/T');
E          = hdf5read(filename,'/arrays/E');
d          = hdf5read(filename,'/arrays/d');
map        = hdf5read(filename,'/arrays/map');
stress     = hdf5read(filename,'/arrays/stress');
visco      = hdf5read(filename,'/arrays/visco');
nT         = params(1);
nE         = params(2);
nd         = params(3);
stress     = cast(stress, 'double');
T          = T - 273.15;

%**** COLOR MAP ****%
%        Dis.         Dif.         GBS          Pei.         MC           Cst. 
mymap = [195,043,037; 251,150,093; 106,193,087; 032,132,199; 101,127,131; 222,138,217]./255;

%**** LOOP ON GRAIN SIZES ****%
for iz=1:nd
    
    %**** GLOBAL INDEX IN THE 3D MAP ****%
    ind     = 1+(iz-1)*nT*nE:nT*nE+(iz-1)*nT*nE;
    map1    = map(ind);
    map2    = reshape(map1,nT,nE)';
    stress1 = stress(ind);
    stress2 = reshape(stress1,nT,nE)';
    visco1  = visco(ind);
    visco2  = reshape(visco1,nT,nE)';
    
    %**** 2D COLORMAP ****%
    figure(iz), clf
    hold on
    imagesc(T,log10(E),map2)
    % Visco contours
    [c,h]=contour(T,log10(E),log10(visco2), [1:50]);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
    clabel(c,h,'FontSize',12,'Color','k')
    set(h, 'Color', 'k', 'LineWidth', 1.0);
    % Stress contours
    [c,h]=contour(T,log10(E),(stress2/1e6), [1 10 100 500 1000]);
    set(h, 'Color', 'w', 'LineWidth', 1.0);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
    clabel(c,h,'FontSize',12,'Color','w')
    colormap(mymap)
    hcb=colorbar; caxis([1 6]), shading flat
    set(gca,'YDir','normal')
    ylabel('Strain rate [s$^{-1}$]','FontSize', 15, 'interpreter', 'latex')
    xlabel('Temperature [$^o$C]'   ,'FontSize', 15, 'interpreter', 'latex')
    title(['Grain size = $10^{', num2str(log10(d(iz)), '%d'), '}$ [m]'],'FontSize', 15 , 'interpreter', 'latex')
    set(gca,'FontSize', 15)
    set(hcb,'YTick',1:6)
    set(hcb,'YTickLabel',{'Dis.','Dif.','GBS','Pei.','MC','Cst.'})
    axis([min(T) max(T) min(log10(E)) max(log10(E))])
    drawnow
    hold off
%     print([path,'./Fig_DefMap', num2str(phase,'%02d'), '_d', num2str(iz,'%02d')], '-dpdf', '-r300')

end