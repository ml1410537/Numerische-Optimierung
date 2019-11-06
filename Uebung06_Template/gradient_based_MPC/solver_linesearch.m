function alpha = solver_linesearch(u, u_prev, dHdu, dHdu_prev, opt)
%SOLVER_LINESEARCH Solve linesearch problem to get stepsize alpha

    % explicit linesearch (see GRAMPC documentation Section 4.3.2 for details)
    delta_u = u - u_prev;
    delta_dHdu = dHdu - dHdu_prev;
    num = delta_u(:)' * delta_u(:);
    den = delta_u(:)' * delta_dHdu(:);
    
    % fallback to LineSearchInit if denominator is close to zero
    if abs(den) < eps
        alpha = opt.LineSearchInit;
        return
    end
    
    alpha = num / den;
    
    % fallback to LineSearchInit if computed stepsize is negative
    if alpha < 0
        alpha = opt.LineSearchInit;
    end
    
    % limit step size to LineSearchMax
    if alpha > opt.LineSearchMax
        alpha = opt.LineSearchMax;
    end
end