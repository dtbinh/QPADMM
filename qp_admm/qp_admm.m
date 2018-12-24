function [x,z,yhat,history] = ...
    qp_admm(P,q,A,l,u,rho,x,z,yhat)
%QP_ADMM -- A QP solver using the ADMM algorithm - James Fleming

VERBOSITY = 1;
MAX_ITERATIONS = 1e5;
T_LIMIT = 30;
EPS_ABS = 1e-6;
EPS_REL = 1e-4;
TAU_INCR = 2;
TAU_DECR = 2;
MU = 10;
ALPH = 1.7;    % relaxation param (1 = normal. [1.5,1.8] = overrelaxation)

% make everything sparse if it isn't already
P = sparse(P);
A = sparse(A);

if nargout == 4
    doHistory = true;
else
    doHistory = false;
end

[m,n] = size(A);

if nargin < 6 || isempty(rho)
    rho = ones(m,1);
end
eqIdx = l == u;
rho(eqIdx) = 1e4*rho(eqIdx);

% From here on, we switch to the preconditioned problem
if false
    % TODO: a decent choice of preconditioner
    [Einv,Ahat,Dinv] = prescale_ruiz(A);
else
    Einv = speye(m);
    Dinv = speye(n);
    Ahat = A;
end

E = inv(Einv);
D = inv(Dinv);
Phat = D'*P*D;
qhat = D'*q;
lhat = E*l;
uhat = E*u;

if nargin ~= 9
    xhat = zeros(n,1);
    zhat = zeros(m,1);
    yhat = zeros(m,1);
%     z = E\zhat;
%     x = D*xhat;
else
    xhat = D\x;
    zhat = E*z;
end

if doHistory
    history = struct();
%     history.x = x;
%     history.z = z;
    history.xhat = xhat;
    history.zhat = zhat;
    history.yhat = yhat;
end

if VERBOSITY >= 1
    fprintf('%3s\t%5s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iteration', 'rho',...
  'r norm', 'eps pri', 's norm', 'eps dual', 'objval');
end

update_factors = true;      % compute factors on first iteration

tStart = tic;

%Atz = Ahat'*zhat;
Aty = Ahat'*yhat;
Atrhoz = Ahat'*(rho.*zhat);

for it=1:MAX_ITERATIONS
    
    % check if Cholesky factors need update
    if update_factors
        [L,~,S] = chol(Phat + Ahat'*(rho.*Ahat), 'lower');
        update_factors = false;
    end
        
    % x update
    %x = (P + rho*(A'*A)) \ (-q - A'*y + rho*A'*z);
%     v = (-qhat - Aty + Ahat'*(rho.*zhat));
%     v = S' * v;
%     v = L \ v;
%     v = L' \ v;
%     xhat = S * v;
   xhat = S * (L' \ (L \ (S' * (-qhat - Aty + Atrhoz))));

    % z update
    %zhat_old = zhat;
    Ax = Ahat*xhat;
    Atrhoz_old = Atrhoz;
    %Ax = ALPH*Ahat*xhat + (1 - ALPH)*zhat;       % overrelaxation
    zhat = min(uhat, max(lhat, (Ax + yhat./rho)));
    Atrhoz = Ahat'*(rho.*zhat);
    
    % y update
    yhat = yhat + rho.*(Ax - zhat);
    Aty = Ahat'*yhat;
    
    % residuals
    %z_old = z;
    %z = E\zhat;
    %x = D*xhat;
    r_norm = norm(Ax - zhat);
    s_norm = norm(Atrhoz_old - Atrhoz);
    
    % new primal dual tolerances - NB: evaluate without preconditioning
    eps_primal = sqrt(m)*EPS_ABS + EPS_REL*max(norm(Ax),norm(zhat));
    eps_dual = sqrt(n)*EPS_ABS + EPS_REL*norm(Aty);
    
    % stopping
    if r_norm < eps_primal && s_norm < eps_dual
        break;
    end
    
    % time limit
    if toc(tStart) > T_LIMIT
        break;
    end
    
    if doHistory
%         history.x = [history.x, x];
%         history.z = [history.z, z];
        history.xhat = [history.xhat, xhat];
        history.zhat = [history.zhat, zhat];
        history.yhat = [history.yhat, yhat];
    end
    
    if VERBOSITY >= 2
        objval = 0.5*xhat'*Phat*xhat + qhat'*xhat;     % objective
        fprintf('%3d\t%4.3f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', it, ...
            rho, r_norm, eps_primal, s_norm, eps_dual, objval);
    end
    
    % update rho
    if r_norm > MU*s_norm
        rho = TAU_INCR*rho;
        update_factors = true;
        Atrhoz = Ahat'*(rho.*zhat);
    elseif s_norm > MU*r_norm
        rho = rho/TAU_DECR;
        update_factors = true;
        Atrhoz = Ahat'*(rho.*zhat);
    end
    
end

if VERBOSITY == 1
    objval = 0.5*xhat'*Phat*xhat + qhat'*xhat;     % objective
    fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', it, ...
        r_norm, eps_primal, s_norm, eps_dual, objval);
end

% Recover real variables
x = D*xhat;
z = E\zhat;
