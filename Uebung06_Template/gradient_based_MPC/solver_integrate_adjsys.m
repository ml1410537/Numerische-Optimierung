function adjvec = solver_integrate_adjsys(tvec, xvec, uvec, p)
%SOLVER_INTEGRATE_ADJSYS Integrate adjoint system in backward time

    % evaluate terminal condition for adjoint state
    adjT = problem_dVdx(xvec(:, end), p);

    % backward time integration of adjoint system
    [~, adj] = ode45(@solver_adjsys, tvec(end:-1:1), adjT, [], tvec, xvec, uvec, p);
    adjvec = adj(end:-1:1, 1:p.Nx)';
end

function out = solver_adjsys(t, adj, tvec, xvec, uvec, p)
%SOLVER_ADJSYS Helper function for integration of adjoint system

    % interpolate state and control
    x = interp1(tvec, xvec', t)';
    u = interp1(tvec, uvec', t)';
    
    % compute adjoint system
    out = -problem_dldx(x, u, t, p) - problem_dfdx_vec(x, u, adj, p);
end