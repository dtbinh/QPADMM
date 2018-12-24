% A test of the qp_admm solver - James Fleming

% load test problem:
load('CVXQP1_M.mat');

% primal problem
[m,n] = size(A);
P_p = Q;
q_p = c;
A_p = [A; speye(n)];
u_p = [ru; ub];
l_p = [rl; lb];

tic;
[x_p,z_p,y_p] = qp_admm(P_p,q_p,A_p,l_p,u_p);
toc;

% tic;
% [x_p,z_p,y_p] = qp_admm_old(P_p,q_p,A_p,l_p,u_p);
% toc;

% % dual problem
% %Q = inv(Q);
% At = [A_p', -A_p'];
% bt = [l_p', -u_p'];
% P_d = At'*(P_p\At);
% q_d = -(bt + q_p'*(P_p\At))';
% A_d = speye(2*size(u_p,1));
% u_d = inf(2*size(u_p,1),1);
% l_d = zeros(2*size(u_p,1),1);
% 
% tic;
% [x_d,z_d,y_d] = qp_admm(P_d,q_d,A_d,l_d,u_d);
% toc;

