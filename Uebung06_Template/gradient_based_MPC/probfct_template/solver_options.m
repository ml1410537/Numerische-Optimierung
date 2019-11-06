function opt = solver_options
%SOLVER_OPTIONS Set options for MPC solver

    % number of discretization steps
    opt.Nhor = 40;
    
    % number of gradient iterations per sampling step
    opt.MaxGradIter = 2;
    
    % fixed step size as fallback for linesearch
    opt.LineSearchInit = 1e-4;
    
    % maximum step size
    opt.LineSearchMax = 0.75;
    
    % shift control trajectory by sampling time
    opt.ShiftControl = true;
end