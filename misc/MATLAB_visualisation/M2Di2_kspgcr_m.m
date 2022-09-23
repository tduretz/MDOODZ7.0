function [y,norm_r,its] = M2Di2_kspgcr_m(A,b,x,L,perm,eps,noisy,SuiteSparse)
% Originally written by Dave A. May - see reference in main text
restart = 20;
max_it  = 1000;
if SuiteSparse==1, r = cs_gaxpy(-A, x, b);
else               r = -A*x + b; end
norm_r  = norm(r);
ncycles = 0;
its     = 0;
rnorm0  = norm_r;
if noisy==3, fprintf(1,'       %1.4d KSP GCR Residual %1.6e %1.6e\n',0,norm_r,norm_r/rnorm0); end
s   = zeros(size(A,1),1);
VV  = zeros(length(x),restart);
SS  = zeros(length(x),restart);
val = zeros(restart  ,1 );
while its<max_it
    for k = 1:restart
        % Preconditionning
        if SuiteSparse==1, s(perm) = cs_ltsolve(L, cs_lsolve(L,r(perm)));
        else               s(perm) = L'\(L\r(perm)); end
        % Action of Jacobian on v
        if SuiteSparse==1, v = cs_gaxpy(A, s, zeros(size(s)));
        else               v = A*s; end
        % Approximation of the Jv product
        for i = 1:k-1
            val(i) = dot(v, VV(:,i));
        end
        for i = 1:k
            v = v - val(i)*VV(:,i);
            s = s - val(i)*SS(:,i);
        end
        r_dot_v = dot(r, v);
        nrm     = sqrt(dot(v, v));
        r_dot_v = r_dot_v/nrm;
        v       = v/nrm;
        s       = s/nrm;
        x       = x + r_dot_v*s;
        r       = r - r_dot_v*v;
        norm_r  = norm(r);
        its     = its+1;
        if noisy==3, fprintf(1,'[%1.4d] %1.4d KSP GCR Residual %1.6e\n',ncycles,its,norm_r/rnorm0);
            fprintf(1,'[%1.4d] %1.4d KSP GCR Residual %1.6e %1.6e\n',ncycles,its,norm_r,norm_r/rnorm0); end
        % Check convergence
        if norm_r < eps*rnorm0, break; end
        VV(:,k) = v;
        SS(:,k) = s;
    end
    if k ~= restart, break; end  % terminated restart loop earlier
    ncycles = ncycles+1;
end
y = x;
end