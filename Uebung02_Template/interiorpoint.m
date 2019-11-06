function sol = interiorpoint(A, b, c, x0, opt)
%INTERIORPOINT Interior point algorithm for linear optimization problems
%
% Vorlesung 'Numerische Optimierung und modellprädiktive Regelung'
%
% Template für Interior Point-Verfahren
%
% Ergänzen Sie die '???'-Einträge
%
% ?bergabeargumente:
%
% A, b, c   Matrizen und Vektoren welche das Optimierungsproblem definieren
% x0        Anfangslösung
% opt       Struktur mit einstellbaren Optionen

    sol = [];
    
    % verify number of arguments
    if nargin < 4
        fprintf('The arguments A, b, c, x0 are mandatory\n');
        return
    end
    if nargin < 5
        opt = [];
    end
    
    % number of optimization variables
    n = size(A, 2);
    % number of equality constraints
    p = size(A, 1);
    
    % verify arguments
    if ~isequal(size(b), [p, 1])
        fprintf('Size of matrix b must be %d x 1\n', p);
        return
    end
    if ~isequal(size(c), [n, 1])
        fprintf('Size of matrix c must be %d x 1\n', n);
        return
    end
    if ~isequal(size(x0), [n, 1])
        fprintf('Size of matrix x0 must be %d x 1\n', n);
        return
    end
    
    % initialize options
    if ~isfield(opt, 'BarrierParamDecreaseFactor')
        opt.BarrierParamDecreaseFactor = 0.1;
    end
    if ~isfield(opt, 'BarrierParamInit')
        opt.BarrierParamInit = 0.1;
    end
    if ~isfield(opt, 'BarrierParamMin')
        opt.BarrierParamMin = 1e-8;
    end
    if ~isfield(opt, 'BoundaryFactor')
        opt.BoundaryFactor = 0.01;
    end
    if ~isfield(opt, 'LineSearchFactor')
        opt.LineSearchFactor = 0.5;
    end
    if ~isfield(opt, 'LineSearchMin')
        opt.LineSearchMin = 1e-8;
    end
    if ~isfield(opt, 'MaxIter')
        opt.MaxIter = 100;
    end
    if ~isfield(opt, 'Tolerance')
        opt.Tolerance = 1e-8;
    end
    
    % initialize state, multipliers and barrier parameter
    sol.x = x0;
    sol.lambda = zeros(p, 1);
    sol.mu = opt.BarrierParamInit * ones(n, 1);
    sol.tau = opt.BarrierParamInit;
    
    % initialize vector of ones
    e = ones(n, 1);
    
    for k = 1 : opt.MaxIter
        sol.iter = k; 
        
        % compute diagonal matrices
        X = diag(sol.x);
        M = diag(sol.mu);

        % determine search direction
        F = [c + A'*sol.lambda - sol.mu; A*sol.x - b; X*M*e - sol.tau*e];
         
        dF = zeros(n+p+n);
        dF(1:n,n+(1:p)) = A';
        dF(n+(1:p),1:n) = A;
        dF(1:n,n+p+(1:n)) = -eye(n);
        dF(n+p+(1:n),1:n) = M;
        dF(n+p+(1:n),n+p+(1:n)) = X;
        
        dF = [zeros(n,n)      A'           -eye(n);
              A           zeros(p,p)    zeros(p,n);
              M           zeros(n,p)            X];
        
        s = - dF \ F;
        
        % extract search directions for x, lambda and mu
        delta_x = s(1:n);
        delta_lambda = s(n+(1:p));
        delta_mu = s(n+p+(1:n));

        % determine step size using backtracking linesearch
        % require that x > bound and mu > bound
        alpha = 1;
        bound = opt.BoundaryFactor * sol.tau;
        feasible = all(sol.x + alpha * delta_x > bound) && all(sol.mu + alpha*delta_mu > bound);
        while ~feasible
            alpha = opt.LineSearchFactor * alpha;
            feasible = all(sol.x + alpha * delta_x > bound) && all(sol.mu + alpha*delta_mu > bound);
            
            if alpha < opt.LineSearchMin && ~feasible
                fprintf('No acceptable step size above %f found\n', opt.LineSearchMin);
                return
            end
        end

        % update variables
        sol.x = sol.x + alpha * delta_x;
        sol.lambda = sol.lambda + alpha * delta_lambda;
        sol.mu = sol.mu + alpha * delta_mu;
        sol.f = c' * sol.x;

        % update barrier parameter
        bar_mu = sol.x'*sol.mu/n;
        sol.tau = max(opt.BarrierParamMin, opt.BarrierParamDecreaseFactor * bar_mu);

        % check convergence
        sol.err = max(F);
        if sol.tau <= opt.BarrierParamMin && sol.err <= opt.Tolerance
            break
        end
    end
end