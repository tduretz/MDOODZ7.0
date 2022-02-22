clear

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
K   = 1.0;
dt = min(dP,dT)^2/K/4.2;
nt  = 5;

SiO2_smooth = SiO2;

for it=1:nt
     SiO2_smooth(2:end-1,2:end-1) = SiO2_smooth(2:end-1,2:end-1) + dt*K*diff(SiO2_smooth(:,2:end-1),2,1)/dT^2 + dt*K*diff(SiO2_smooth(2:end-1,:),2,2)/dP^2;
end

[~,T600] = min( abs(T-273-600) );

fileID  = fopen(['SiO2_nsm', num2str(nt, '%03d'), '.dat'],'w');
fwrite(fileID, SiO2_smooth, 'double');
fclose(fileID);

figure(1)
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
subplot(3,1,3), title('Density profiles')
plot(P/1e9,SiO2(T600,:),'k', P/1e9,SiO2_smooth(T600,:),'b')
xlabel('P [GPa]'), ylabel('rho')


