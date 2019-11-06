function sol = simplex(A, b, c, Basis)
%SIMPLEX Simplex-algorithm for linear optimization problems
%
% Vorlesung 'Numerische Optimierung und modellprädiktive Regelung'
%
% Template für Simplex-Verfahren
%
% Ergänzen Sie die '???'-Einträge
%
% ?bergabeargumente:
%
% A, b, c   Matrizen und Vektoren welche das Optimierungsproblem definieren
% Basis     Vektor welcher die Basisvariablen definiert

    sol = [];
    
    % verify number of arguments
    if nargin < 3
        fprintf('The arguments A, b, c are mandatory\n');
        return
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

    % initialization phase
    if nargin < 4
        E = diag((b >= 0) - (b < 0));
        solinit = simplex([A E], b, [zeros(n,1); ones(p,1)], n+(1:p)); 
        if ~all(solinit.x(n+(1:p)) == 0)
            fprintf('Infeasible optimization problem\n');
            return
        end
        Basis = solinit.Basis;
    end    
    
    % determine nonbasis from basis
    Nonbasis = setdiff(1:n, Basis);   
    
    k = 0;
    while 1
        k = k + 1;
        
        % separate matrix A in B and N
        B = A(:,Basis);
        N = A(:,Nonbasis);

        % separate vector c in cB and cN
        cB = c(Basis);
        cN = c(Nonbasis);

        % compute xB
        xB = B\b;
        
        if ~all(xB >= 0)
            fprintf('Invalid basis\n');
            return
        end

        % compute multipliers for equalities
        lambda = -B'\cB;
        
        % compute multipliers for active inequalities
        muN = cN + N'*lambda;
        
        % optimal solution found?
        if all(muN >= 0)
            if nargin < 4
                fprintf('Optimal solution found\n');
            end
            % set optimal basis and nonbasis
            sol.Basis = Basis;
            % set optimal x
            sol.x = zeros(n, 1);
            sol.x(Basis) = xB;
            % set optimal lambda
            sol.lambda = lambda;
            % set optimal mu
            sol.mu = zeros(n, 1);
            sol.mu(Nonbasis) = muN;
            % set optimal cost
            sol.f = c' * sol.x;
            % set number of iterations
            sol.iter = k;
            return
        end

        % determine index of minimal multiplier
        [~, si] = min(muN);
        s = Nonbasis(si);
        % compute direction d
        d = B\A(:,s);
        % is problem unconstrained in the direction of d?
        if all(d <= 0)
            fprintf('Unconstrained problem\n');
            return
        end

        % determine index of variable that reaches the constraint first
        xsp = inf;
        for i = 1 : p
            if (d(i)>0) && (xB(i)/d(i)<xsp)
                xsp = xB(i)/d(i);
                ri = i;
            end
        end
        r = Basis(ri);

        % switch basis and nonbasis variables
        Basis(ri) = s;
        Nonbasis(si) = r;
    end
end
