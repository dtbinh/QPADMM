% A test of the qp_admm solver - James Fleming

% load test problem:
load('CVXQP1_M.mat');

Aqp = [A;-A];
bqp = [ru; -rl];

% tic;
% [x,fval,exitflag,output,lambda] = quadprog(Q,c,Aqp,bqp,[],[],lb,ub);
% toc;

[m,n] = size(A);
P = Q;
q = c;
Aadmm = [A; speye(n)];
u = [ru; ub];
l = [rl; lb];

% tic;
% [x_p,z_p,y_p] = qp_admm_noeq(P,q,Aadmm,l,u);
% toc;

tic;
[x_p,z_p,y_p,mu_p] = qp_admm(P,q,Aadmm,l,u);
toc;

% % perturbed initial guess
% x_g = x_admm + 1e-6*randn(size(x_admm));
% z_g = z_admm + 1e-6*randn(size(z_admm));
% y_g = y_admm + 1e-6*randn(size(y_admm));
% 
% tic;
% [x_admm2,z_admm2,y_admm2] = qp_admm(P,q,Aadmm,l,u,1e3,x_g,z_g,y_g);
% toc;