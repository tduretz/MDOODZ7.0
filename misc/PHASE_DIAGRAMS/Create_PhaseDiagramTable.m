clear, clc, figure(1), clf, colormap(jet) % Plotting densities, th. exp & compressibility coeff. calculated with PERPLE_X
fprintf('# =========================================== #\n');
fprintf('#         CREATE PHASE DIAGRAM TABLES         #\n')
fprintf('# =========================================== #\n');
%% Miscellaneous
nT = 1500;                                                                                                  % Resolution for temperature []
nP = 1500;                                                                                                  % Resolution for pressure    []
fname = '/Users/lcandiot/Software/mdoodz4.5_repo/SOURCE/PHASE_DIAGRAMS/HawaiianPyrolite_HR_rho_bin.dat';   % Path to local copy of mdoodz repository
write = 1;                                                                                                  % Writing binary
%% Read from table
fprintf('Loading data ...\n');
tic
load ~/Dropbox/Vangelis_Lorenzo_Jenadi/work/Jenadi_stx_HR_1.tab;            % Perplex table
fprintf('Time for loading data: %.2f s\n',toc);
Dat = Jenadi_stx_HR_1;
%% Initialisation
T   = Dat(:,1); T2d   = reshape(T  ,nT,nP)';         % Massage arrays
P   = Dat(:,2); P2d   = reshape(P  ,nT,nP)';
rho = Dat(:,3); 
alp = Dat(:,4); 
bet = Dat(:,5); 
%% Check for NaNs in database
fprintf('Checking for NaNs in dataset ...\n');
tic
rho = fillmissing(rho,'linear');
alp = fillmissing(alp,'linear');
bet = fillmissing(bet,'linear');
fprintf('Time for NaN check: %.2f s\n', toc);
% Reshape
rho2d = reshape(rho,nP,nT)'; %rho2d = flipud(rho2d);
bet2d = reshape(bet,nT,nP);
alp2d = reshape(alp,nT,nP);
T_plot = unique(T); P_plot = unique(P);
%% Visualisation
figure(1),clf
imagesc(T_plot,P_plot,rho2d), colorbar
axis xy
hold on
[c,hc] = contour(T_plot-273.15,P_plot/1e8,rho2d,[3000:20:4000]);
set(hc,'Color','k')
fprintf('Please check the data field before I will write to disk!\n');
pause
%% Write to binary
if(write==1), fprintf('Writing density ...\n'); tic; fid=fopen( fname,'w' ); fwrite(fid, rho,'double',0,'n'); fclose(fid); fprintf('Time to write densities: %.2f s\n',toc); end