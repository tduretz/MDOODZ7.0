clear

step = 1;

% Data from MDOODZ:
path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE/';
filename = [path 'Stokes_04cpu_step', num2str(step,'%02d'), '_iter22.gzip.h5'];
[ MA, MB, MC, MD, rSt, rPt, fSt, fPt ] = ExtractDecoupledMatrixFromMDoodz( filename );

% Data from MDOODZ:
path = '/Users/imac/REPO_GIT/MDOODZ6.0/SOURCE/';
filename = [path 'Jacobian_04cpu_step', num2str(step,'%02d'), '_iter22.gzip.h5'];
[ JA, JB, JC, JD, rJSt, rJPt, fJSt, fJPt ] = ExtractDecoupledMatrixFromMDoodz( filename );

D = 0*diag(size(MC,1)-1,size(MC,1)-1);

% Matrices
K = [MA MB; MC D ] ;
J = [JA JB; JC D ] ;


diff = K-J;
[i,j,s] = find(diff);
s(abs(s)<1e-6) = 0;
diff = sparse(i,j,s);
min(diff(:))
max(diff(:))
% diff(diff<1e-4) = 0;

figure(1)
subplot(131), spy(K)
subplot(132), spy(J)
subplot(133), spy(diff)

% fprintf('Picard')
% K(1,:)
% fprintf('Jacobian')
% J(1,:)
% 
% fprintf('Matlab')
% M1(699,:)
% fprintf('M2Di2')
% M2(699,:)
% 
% fprintf('Matlab')
% M1(2600,:)
% fprintf('M2Di2')
% M2(2600,:)

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