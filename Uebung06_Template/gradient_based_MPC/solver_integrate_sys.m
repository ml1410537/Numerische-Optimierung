function [xvec, J] = solver_integrate_sys(x0, tvec, uvec, p)
%SOLVER_INTEGRATE_SYS Integrate system in forward time and compute cost

        % integrate system and cost
        [~, x] = ode45(@solver_sys, tvec, [x0; 0], [], tvec, uvec, p);
        xvec = x(:, 1:p.Nx)';
        
        % add terminal cost
        J = x(end, end) + problem_Vfct(xvec(end, 1:p.Nx)', p);
end

function out = solver_sys(t, x, tvec, uvec, p)
%SOLVER_SYS Helper function for integration of system and cost

    % interpolate control
    u = interp1(tvec, uvec', t)';
    
    % compute system and cost
    out = [problem_ffct(x(1:p.Nx), u, p);
           problem_lfct(x(1:p.Nx), u, t, p)];
end