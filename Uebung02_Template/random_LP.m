function [A, b, c, x0] = random_LP(n, p)
%RANDOM_LP Generate random linear program with n variables, p constraints
    % cf. Boyd, "Convex Optimization", p. 575
    %
    % Vorlesung 'Numerische Optimierung und modellprädiktive Regelung'
    %
    if nargin < 2
        % We take twice as many variables as constraints
        p = round(n / 2);
    end
    % The elements of A are iid N(0,1)
    A = randn(p, n);
    % The elements of x0 are iid U(0,1)
    x0 = rand(n, 1);
    % This ensures primal feasibility
    b = A * x0;
    % We compute a vector iid N(0,1)
    z = randn(p, 1);
    % and a vector iid U(0,1)
    s = rand(n, 1);
    % and then take
    c = A' * z + s;
    % which ensures dual feasibility
end

