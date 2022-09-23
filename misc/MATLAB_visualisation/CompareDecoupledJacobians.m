clear

Mc = load('JacobianOK');
Md = load('Jacobian');

% Matrices
M1 = [Mc.MA Mc.MB; Mc.MC 0*Mc.MD ];
M2 = [Md.MA Md.MB; Md.MC 0*Md.MD ];

% RHSs
% R1 = Mc.L;
% R2 = [Md.rSt; Md.rPt];

diff = M1-M2;

figure(1)
subplot(131), spy(M1)
subplot(132), spy(M2)
subplot(133), spy(diff)

diff(1,:)

% norm(R1-R2)

% Mc = load('CoupledMatrix_BENALLAL');
% Md = load('CoupledMatrix_BENALLAL2');
% 
% % Matrices
% M1 = Mc.M;
% M2 = Md.M;
% 
% % RHSs
% R1 = Mc.L;
% R2 = Md.L;
% 
% figure(1)
% subplot(131), spy(M1)
% subplot(132), spy(M2)
% 
% 
% Mdiff = M1-M2;
% Mdiff = 0;
% subplot(133), spy(Mdiff)
% 
% norm(R1-R2)