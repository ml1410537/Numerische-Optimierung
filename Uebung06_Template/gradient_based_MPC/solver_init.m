function [vec, sol] = solver_init(p, opt, Tsim)
%SOLVER_INIT Initializes structures for MPC solver

    % vector of MPC solution
    vec.t = 0 : p.dt : Tsim;
    Nsim = length(vec.t);
    vec.x = nan * ones(p.Nx, Nsim);
    vec.u = nan * ones(p.Nu, Nsim);
    vec.adj = nan * ones(p.Nx, Nsim);
    vec.J = nan * ones(1, Nsim);
    
    % solution of OCP
    sol.t = linspace(0, p.Thor, opt.Nhor);
    sol.x = repmat(p.x0, 1, opt.Nhor);
    sol.u = repmat(p.u0, 1, opt.Nhor);
    
    sol.u_prev = sol.u;
    sol.dHdu = 0 * sol.u;
    sol.dHdu_prev = 0 * sol.u;
end

