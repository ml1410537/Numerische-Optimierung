function [G, c, Aeq, beq, C, d] = lqr_problem(A, B, Q, R, S, x0, N, xmin, xmax, umin, umax)
%LQR_PROBLEM Construct QP for the linear-quadratic problem
%
% min_u  \sum_{k=0}^{N-1} (x_k' Q x_k + u_k' R u_k) + x_N' S x_N
% s.t.   x_{k+1} = A x_k + B u_k
%        xmin <= x_k <= xmax
%        umin <= u_k <= umax
%
% in the form of
%
% min_{z}  0.5 z' G z + c' z
% s.t.     Aeq z = beq
%          C z <= d

    Nx = size(A, 2);
    Nu = size(B, 2);
    
    % cost function
    G = zeros((Nx+Nu)*N);
    c = zeros((Nx+Nu)*N, 1);
    for i = 1 : N
        if i ~= N
            G((Nx+Nu)*(i-1)+1:(Nx+Nu)*i-1,(Nx+Nu)*(i-1)+1:(Nx+Nu)*i-1) = Q;
            G((Nx+Nu)*i,(Nx+Nu)*i) = R;
        else
            G((Nx+Nu)*(i-1)+1:(Nx+Nu)*i-1,(Nx+Nu)*(i-1)+1:(Nx+Nu)*i-1) = S;
            G((Nx+Nu)*i,(Nx+Nu)*i) = R;
        end
    end
    
    % equality constraints
    Aeq = zeros(Nx*N, (Nx+Nu)*N);
    beq = zeros(Nx*N, 1);
    for i = 1 : N
       if i == 1
           Aeq(Nx*(i-1)+1:Nx*i,(Nx+Nu)*(i-1)+1:(Nx+Nu)*i) = [eye(Nx) -B];
           beq(1:Nx) = [A*x0];
       else
           Aeq(Nx*(i-1)+1:Nx*i,(Nx+Nu)*(i-2)+1:(Nx+Nu)*i) = [-A zeros(Nx,Nu) eye(Nx) -B];
       end
    end
    
    % add inequality constraints if xmin, xmax, umin, umax are given
    if nargin == 11
        C = zeros(2*(Nx+Nu)*N, (Nx+Nu)*N);
        d = zeros(2*(Nx+Nu)*N, 1);
        for i = 1 : N
            C(2*(Nx+Nu)*(i-1)+1:2*(Nx+Nu)*i,(Nx+Nu)*(i-1)+1:(Nx+Nu)*i) = [eye(Nx) zeros(Nx,Nu);-eye(Nx) zeros(Nx,Nu);zeros(Nx) [1;-1]];
            d(2*(Nx+Nu)*(i-1)+1:2*(Nx+Nu)*i) = [xmax;-xmin;umax;-umin];
        end
    else
        C = [];
        d = [];
    end

end

