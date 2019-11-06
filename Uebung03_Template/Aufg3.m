% ------ Initialisierung ------ 
T_hor = 3;
N = 30;
delta_t =T_hor/N;


S = diag([1,0]);
Q = diag([1,0])*delta_t;
R = 0.01*delta_t;
x0 = [1;0];
u = 0;

A = [1 delta_t;0 1];
B = [0;delta_t];
Nx = size((A),1);
Nu = size((B),2);


% ------ Aufg 3.2 Riccati ------ 
disp("================== Aufg 3.2 Riccati ================== ")
[G, c, Aeq, beq] = lqr_problem(A, B, Q, R, S, x0, N);
[x_ri, u, iter] = riccati(A, B, Q, R, S, x0, N);
[x_lqr,~,~,output] = quadprog(G, c, Aeq, beq, [], [],[x0;0;zeros(length(c)-3,1)]);
%sol_lqr = active_set_qp(G, c, Aeq, beq);
%h(2) = stairs(x_lqr);
%x1(1) = x0(1);
%for i=2:N
%   x1(i) = x_lqr((3*(i-1)+1),1); 
%end

% ------ Aufg 3.3 Aktive-Restriktionen-Verfahren ------ 
xmax = [1;1];
xmin = [-1;-1];
umax = 1;
umin =-1;
[G, c, Aeq, beq, C, d] = lqr_problem(A, B, Q, R, S, x0, N, xmin, xmax, umin, umax);
disp("================== Aufg 3.3 Aktive-Restriktionen-Verfahren ================== ")
sol.x  = active_set_qp(G, c, Aeq, beq, C, d);

% ------ Aufg 3.4 Interiorpoint_qp ------ 
disp("================== Aufg 3.4 Interiorpoint_qp ================== ")
opt.MaxIter = 200;
opt.BarrierParamDecreaseFactor = 0.001;
sol = interiorpoint_qp(G, c, Aeq, beq, 0.1*ones(size(c)), opt);
reshape(sol.x, Nx+Nu, N)
% ------ Aufg 3.5 Polygon ------ 
disp("================== Aufg 3.5 Polygon ================== ")
solve_polygon(6)