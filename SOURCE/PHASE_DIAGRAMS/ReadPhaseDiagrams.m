clear
% Check SiO2
fileID = fopen('SiO2.dat','r');
nT     = 1300;
nP     = 1000;           
Tmin    = 298+1e-3;            % K
Tmax    = 1298+300;
Pmin    = 0.049/10*1e9;        % Pa
Pmax    = 99.95/10*1e9;
dT      = (Tmax-Tmin)/(nT-1);
dP      = (Pmax-Pmin)/(nP-1);
T       = Tmin:dT:Tmax;
P       = Pmin:dP:Pmax;
SiO2   = fread(fileID, nT*nP, 'double');
fclose(fileID);
SiO2   = reshape(SiO2, nT, nP);

figure(1), clf
subplot(1,2,1), title(SiO2)
imagesc(T-273, P/1e9, SiO2'), axis xy, title('Si02')
xlim([400 1300])
ylim([0 5.5])
colorbar

% Check Rhyolite
fileID = fopen('Rhyolite.dat','r');  
nT      = 1249;
nP      = 1249;                 % kbar
Tmin    = 373.0;               % K
Tmax    = 1373.0;
Pmin    = 10.13e6;             % Pa
Pmax    = 5.5816e9;
dT      = (Tmax-Tmin)/(nT-1);
dP      = (Pmax-Pmin)/(nP-1);
T       = Tmin:dT:Tmax;
P       = Pmin:dP:Pmax;
Rhyo    = fread(fileID, nT*nP,'double');
fclose(fileID);
Rhyo    = reshape(Rhyo, nT, nP);

figure(1),
subplot(1,2,2)
imagesc(T, P/1e9, Rhyo'), axis xy, title('Rhyolite')
xlim([400 1300])
ylim([0 5.5])
colorbar


