clear

%% For comparison, load 1D @ 600C
fileID  = fopen('QuartzCoesite600C.bin', 'r');
nP      = 10000;
Pmin    = 1e5 /1e9*10; % kbar
Pmax    = 1e10/1e9*10;
dP      = (Pmax-Pmin)/(nP-1);
P1D     = Pmin:dP:Pmax;
SiO21D  = fread(fileID, nP, 'double');
fclose(fileID);

fileID  = fopen('QuartzCoesite600C_smoothed.bin', 'r');
nP      = 10000;
Pmin    = 1e5 /1e9*10; % kbar
Pmax    = 1e10/1e9*10;
dP      = (Pmax-Pmin)/(nP-1);
P1D     = Pmin:dP:Pmax;
SiO21Dsm= fread(fileID, nP, 'double');
rhoc    = 1.0000000000e+43;
fclose(fileID);

% Check SiO2
fileID  = fopen('SiO2.dat','r');
nT      = 675;
nP      = 2500;           
Tmin    = 398+1e-3;            % K
Tmax    = 800+273;
Pmin    = 0.0090;        % Pa
Pmax    = 49.9890; %/10*1e9
dT      = (Tmax-Tmin)/(nT-1);
dP      = (Pmax-Pmin)/(nP-1);
T       = Tmin:dT:Tmax;
P       = Pmin:dP:Pmax;
SiO2    = fread(fileID, nT*nP, 'double');
fclose(fileID);
SiO2    = reshape(SiO2, nT, nP);

% Diffusion smoothing
K           = 1.0;
dt          = min(dP,dT)^2/K/4.2;
nt          = 100;
SiO2_smooth = SiO2;
for it=1:nt
     SiO2_smooth(2:end-1,2:end-1) = SiO2_smooth(2:end-1,2:end-1) + dt*K*diff(SiO2_smooth(:,2:end-1),2,1)/dT^2 + dt*K*diff(SiO2_smooth(2:end-1,:),2,2)/dP^2;
end

[~,T600] = min( abs(T-273-600) );

fileID  = fopen(['SiO2_nsm', num2str(nt, '%03d'), '.dat'],'w');
fwrite(fileID, SiO2_smooth, 'double');
fclose(fileID);

figure(1), clf
subplot(3,1,1), title(SiO2)
imagesc(T-273, P/1e9, SiO2'), axis xy, title('Si02')
hold on, plot((600)*ones(size(P)), P/1e9, 'w', 'LineWidth', 2)
% xlim([400 1300])
% ylim([0 5.5])
colorbar
subplot(3,1,2), title(SiO2_smooth)
imagesc(T-273, P/1e9, SiO2_smooth'), axis xy, title('Si02')
hold on, plot((600)*ones(size(P)), P/1e9, 'w', 'LineWidth', 2)
% xlim([400 1300])
% ylim([0 5.5])
colorbar
subplot(3,1,3), title('Density profiles'), hold on
plot(P,SiO2(T600,:),'k', P,SiO2_smooth(T600,:),'b')
plot(P1D,SiO21D,'-.r')
plot(P1D,SiO21Dsm*rhoc,'r')
xlim([25 30])
xlabel('P [kbar]'), ylabel('rho')


