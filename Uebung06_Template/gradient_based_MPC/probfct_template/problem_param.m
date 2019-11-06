function p = problem_param
%PROBLEM_PARAM Set parameters for MPC problem
    
    % number of states and control
    p.Nx = 10;
    p.Nu = 3;
    
    % initial values
    p.x0 = [0;0;0;0;0;0;0;0;0;0];
    p.u0 = [0;0;0];
    
    % desired values
    p.xdes = [1;0;1;0;0;0;0;0;0;0];
    p.udes = [1;0;0];
    
    % weighting matrices
    p.S = diag([1,0.1,1,0.1,1,0.1,1,0.1,1,0.1]);
    p.Q = diag([1,0.1,1,0.1,1,0.1,1,0.1,1,0.1]);
    p.R = diag([0.01,0.01,0.01]);

    % control limits
    p.umin = [-2;-2;-2];
    p.umax = [2;2;2];
    
    % time horizon
    p.Thor = 1.5;
    % sampling time
    p.dt = 2e-3;
end
