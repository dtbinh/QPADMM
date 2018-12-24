function [x,z,y,history] = qp_admm(P,q,A,l,u,rho,x,z,y)
%QP_ADMM - A QP solver using the ADMM algorithm - James Fleming

VERBOSITY = 1;
MAX_ITERATIONS = 1e5;
T_LIMIT = 30;
EPS_ABS = 1e-9;
EPS_REL = 1e-6;
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
nEq = sum(eqIdx);
rho(eqIdx) = 1e4*rho(eqIdx);

if nargin ~= 9
    x = zeros(n,1);
    z = zeros(m,1);
    y = zeros(m,1);
%     % start from eq constrained QP solution?
%     Aeq = A(eqIdx,:);
%     beq = u(eqIdx,:);
%     res = [P,Aeq';Aeq,zeros(nEq)] \ [-q; beq];
%     x = res(1:n);
%     z = A*x;
%     y = zeros(m,1);
end

if doHistory
    history = struct();
    history.x = x;
    history.z = z;
    history.y = y;
end

if VERBOSITY >= 1
    fprintf('%3s\t%5s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iteration', 'rho',...
  'r norm', 'eps pri', 's norm', 'eps dual', 'objval');
end

update_factors = true;      % compute factors on first iteration

tStart = tic;

Atz = A'*z;
Aty = A'*y;

for it=1:MAX_ITERATIONS
    
    % check if Cholesky factors need update
    if update_factors
        [L,~,S] = chol(P + (A'*sparse(diag(rho))*A), 'lower');
        update_factors = false;
    end
        
    % x update
    %x = (P + rho*(A'*A)) \ (-q - A'*y + rho*A'*z);
    v = (-q - Aty + A'*(rho.*z));
    v = S' * v;
    v = L \ v;
    v = L' \ v;
    x = S * v;

    % z update
    z_old = z;
    %Ax = A*x;
    Ax = ALPH*A*x + (1 - ALPH)*z;       % relaxation
    z = min(u, max(l, (Ax + y./rho)));
    
    % y update
    y = y + rho.*(Ax - z);
    Aty = A'*y;
    
    % residuals
    r_norm = norm(Ax - z);
    s_norm = norm(A'*(rho.*z_old) - A'*(rho.*z));
    
    % new primal dual tolerances
    eps_primal = sqrt(m)*EPS_ABS + EPS_REL*max(norm(Ax),norm(z));
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
        history.x = [history.x, x];
        history.z = [history.z, z];
        history.y = [history.y, y];
    end
    
    if VERBOSITY >= 2
        objval = 0.5*x'*P*x + q'*x;     % objective
        fprintf('%3d\t%4.3f\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', it, ...
            rho, r_norm, eps_primal, s_norm, eps_dual, objval);
    end
    
    % update rho
    if r_norm > MU*s_norm
        rho = TAU_INCR*rho;
        update_factors = true;
    elseif s_norm > MU*r_norm
        rho = rho/TAU_DECR;
        update_factors = true;
    end
    
end

if VERBOSITY == 1
    objval = 0.5*x'*P*x + q'*x;     % objective
    fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', it, ...
        r_norm, eps_primal, s_norm, eps_dual, objval);
end
