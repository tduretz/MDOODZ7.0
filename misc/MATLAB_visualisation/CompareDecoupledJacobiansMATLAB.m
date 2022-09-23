clear

Mc = load('Jacobian_MATLAB');
Md = load('Jacobian');

% Matrices
M1 = [Mc.J1 Mc.grad1; Mc.div1 0*Md.MD ];
M2 = [Md.MA Md.MB; Md.MC 0*Md.MD ];

% Kill dirichlet row/cols from the M2Di code
west = Mc.NumVx(1,:);
east = Mc.NumVx(end,:);
south = Mc.NumVyG(:,1);
north = Mc.NumVyG(:,end);
all = [west(:); east(:); south(:); north(:)];
M1(all,:) = [];
M1(:,all) = [];

% RHSs
% R1 = Mc.L;
% R2 = [Md.rSt; Md.rPt];

diff = abs(M1-M2);
diff(diff<1e-4) = 0;

figure(1)
subplot(131), spy(M1)
subplot(132), spy(M2)
subplot(133), spy(diff)

M1(1,:)
M2(1,:)

M1(2600,:)
M2(2600,:)

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