% A test of the qp_admm solver - James Fleming

% load test problem:
load('QSCSD6.mat');

Aqp = [A;-A];
bqp = [ru; -rl];

tic;
[x,fval,exitflag,output,lambda] = quadprog(Q,c,Aqp,bqp,[],[],lb,ub);
toc;

% Primal problem
[m,n] = size(A);
P = Q;
q = c;
An = [A; speye(n)];
u = [ru; ub];
l = [rl; lb];

tic;
[x_admm,z_admm,y_admm] = qp_admm(P,q,An,l,u);
toc;

% Dual problem
