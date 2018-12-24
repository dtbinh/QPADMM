function [Einv,Ahat,Dinv] = prescale_new(A)
%PRESCALE_RUIZ Ruiz prescaling for preconditioning of matrix A
% See "A scaling algorithm to equilibrate both rows and columns norms
% in matrices" by Daniel Ruiz

MAX_ITER = 100;
EPS = 1e-3;

[m,n] = size(A);

Einv = speye(m);
Dinv = speye(n);

Ahat = A;

for k=1:MAX_ITER
    
    %absA = abs(Ahat);
    dr = sqrt(sqrt(sum(Ahat.^2, 2)));
    dc = sqrt(sqrt(sum(Ahat.^2, 1)));
    
    if abs(1-max(dr)) < EPS && abs(1-max(dc)) < EPS
        %break;
    end
    
    Ahat = dr.\(Ahat./dc);
    Einv = Einv./dr;
    Dinv = Dinv./dc;    
    
end

end

