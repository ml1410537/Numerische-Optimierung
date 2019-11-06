function vec = solver_main(figNr, Tsim)
%SOLVER_MAIN Model predictive control using projected gradient method
%
% Vorlesung "Numerische Optimierung und modellprädiktive Regelung"
% Sommersemester 2019
% Template für MPC mit projiziertem Gradientenverfahren
    % load parameters
    p = problem_param;
    % load options
    opt = solver_options;
    
    if nargin < 1
        figNr = 1;
    end    
    if nargin < 2
        Tsim = p.Thor;
    end
    
    % init structures
    [vec, sol] = solver_init(p, opt, Tsim);
    x0 = p.x0;
    
    % MPC loop
    for k = 1 : length(vec.t)        
        % gradient loop
        for i = 1 : opt.MaxGradIter        
           % TODO implement projected gradient method
           % compute sol.x, sol.adj, sol.u, ...
           sol.adj =  solver_integrate_adjsys(sol.t, sol.x, sol.u, p);
           alpha = solver_linesearch(sol.u,sol.u_prev,sol.dHdu,sol.dHdu_prev,opt);
           sol.dHdu_prev = sol.dHdu;
           sol.dHdu = solver_gradient(sol.t, sol.x, sol.u, sol.adj, p);
           sol.u_prev = sol.u;
           sol.u = solver_project(sol.u_prev - alpha * sol.dHdu, p); 
        end
        
        % forward time integration of system and cost
        [sol.x, sol.J] = solver_integrate_sys(x0, sol.t, sol.u, p);
        
        % store results
        vec.x(:,k) = x0;
        vec.u(:,k) = sol.u(:,1);
        vec.adj(:,k) = sol.adj(:,1);
        vec.J(k) = sol.J;
        
        % go to next sampling step
        x0 = interp1(sol.t, sol.x', p.dt)';
               
        % shift control by sampling time
        if opt.ShiftControl
            sol.u = solver_shiftcontrol(sol.t, sol.u, p.dt);
            sol.uprev = solver_shiftcontrol(sol.t, sol.u_prev, p.dt);
            sol.dHdu = solver_shiftcontrol(sol.t, sol.dHdu, p.dt);
            sol.dHdu_prev = solver_shiftcontrol(sol.t, sol.dHdu_prev, p.dt);
        end
        
        % plot results every 10 steps
        if mod(k, 10) == 0
            figure(figNr);
            solver_plot(vec, sol, vec.t(k));
        end
    end
end
