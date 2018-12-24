% A toy QP problem for testing a QP algorithm - James Fleming
badparam = 1e1; % ADMM struggles as this gets large (>20)

H = [2 -1; -1 5]; 
f = [5; 2];
A = [badparam 1; -1 0];
b = [1; 1];

% NB: The problem seems to be that cond(A'*A) becomes large

P = H;
q = f;
u = b;
l = -inf(size(b));

tic;
[x,fval,exitflag,output,lambda] = quadprog(H,f,A,b);
toc;

% solve preconditioned problem
tic;
[x_p,z_p,y_p,history] = qp_admm(P,q,A,l,u);
toc;

% some plots to see what's going on
figure(101); hold on; cla;

x = linspace(min(history.x(1,:)), max(history.x(1,:)));
y = linspace(min(history.x(2,:)), max(history.x(2,:)));
[X,Y] = meshgrid(x,y);
Z = 0.5*(X.*H(1,1).*X + Y.*H(2,1).*X + X.*H(1,2).*Y + Y.*H(2,2).*Y) ...
    + f(1).*X + f(2).*Y;
feasible = (A(1,1)*X + A(1,2)*Y <= b(1)) & (A(2,1)*X + A(2,2)*Y <= b(2));
Z(~feasible) = nan;
contourf(X,Y,Z,10);

plot(history.x(1,:), history.x(2,:), 'kx');
plot(history.x(1,:), history.x(2,:), 'r--');