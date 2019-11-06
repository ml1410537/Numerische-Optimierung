function sol = active_set_qp(G, c, A, b, C, d)
%ACTIVE_SET Active set method for solving quadratic programs
% min  0.5*x'*Q*x + c'*x
% s.t. A*x = b
%      C*x <= d
%
%
% Vorlesung 'Numerische Optimierung und modellprädiktive Regelung'
%
% Template für Interior Point-Verfahren
%
% Ergänzen Sie die '???'-Einträge
%
% ?bergabeargumente:
%
% G,c  Definition der quadratischen Kostenfunktion
% A,b  Definition der linearen Gleichungsbeschränkungen
% C,d  Definition der linearen Ungleichungsbeschränkungen

    if nargin < 3
        A = [];
        b = [];
    end
    if nargin < 5
        C = [];
        d = [];
    end
    
    % tolerance for optimality
    tolerance = 1e-12;

    % initialize solution
    sol.iter = 0;
    
    % dimension of optimization variable x
    n = length(c);
    % number of equality constraints
    p = length(b);
    % number of inequality constraints
    q = length(d);

    if q == 0
        if p == 0
            % solve unconstrained problem
            sol.x = quadprog(G ,c);
            sol.lambda = [];
        else
            % solve equality-constrained problem
            y = [G A';A zeros(p)]\[-c;b];        
            sol.x = y(1:n);
            sol.lambda = y(n+1:n+p);
        end
        return
    end
    
    % inequality-constrained problem
    % active set method
    
    % initial value
    xk = zeros(n,1);
    % initial active constraint set
    active = find(C * xk >= d);
    % initial inactive constraint set
    inactive = setdiff(1:q, active);
    
    while 1
        % count number of iterations
        sol.iter = sol.iter + 1;
        
        % define equality constraints
        if p > 0
            A_act = [A;C(active,:)];
            b_act = [b - A*xk; d(active) - C(active,:) * xk];
            %b_act = [zeros(p+length(active),1)];
        else
            A_act = [C(active,:)'];
            b_act = [d(active) - C(active,:) * xk];
        end
        % solve subproblem by calling active_set_qp 

        sol_act = active_set_qp(G, (G*xk+c), A_act, b_act);
        sk = sol_act.x;
        mu = sol_act.lambda(p + (1 : length(active)));
        
        if all(abs(sk) < tolerance)
            % minimum found
            if all(mu >= 0)
                % optimal solution found
                fprintf('Optimal solution found\n');
                sol.x = xk;
                sol.lambda = sol_act.lambda(1:p);
                sol.mu = mu;
                return
            else
                % find constraint with minimal multiplier
                [~, idx] = min(mu);
                i = active(idx);
                % deactivate active constraint
                fprintf('Deactivate constraint %i\n', i);
                active = setdiff(active, i);
                inactive = union(inactive, i);
            end
        else
            % determine step length
            alpha = 1;
            i = 0;
            for j = inactive
                if C(j, :) * sk > 0
                    % determine step length for inactive constraint j
                    alpha_j = (d(j)-C(j, :)*xk)/(C(j, :)*sk);
                    if alpha_j < alpha
                        alpha = alpha_j;
                        i = j;
                    end
                end
            end
            if alpha < 1
                % activate constraint
                fprintf('Activate constraint %i\n', i);
                active = union(active, i);
                inactive = setdiff(inactive, i);
            end
            
            % go to next step
            xk = xk + alpha * sk;        
        end
    end
end
