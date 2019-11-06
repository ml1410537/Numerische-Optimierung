function [x, u, k] = riccati(A, B, Q, R, S, x0, N, umin, umax)
%RICCATI Solution of the linear-quadratic problem
%
% min_u  \sum_{k=0}^{N-1} (x_k' Q x_k + u_k' R u_k) + x_N' S x_N
% s.t.   x_{k+1} = A x_k + B u_k
%
% using the discrete-time Riccati equation

    Nx = size(A, 2);
    Nu = size(B, 2);

    % compute Riccati matrix P_k backwards
    P = zeros(Nx, Nx, N+1);
    P(:, :, N+1) = S;
    for k = N : -1 : 1
        Pk = P(:, :, k+1);
        P(:, :, k) = (Q+A'*Pk*A)-(B'*Pk*A)'*inv(R+B'*Pk*B)*(B'*Pk*A);
    end

    % compute control matrix K_k
    K = zeros(Nu, Nx, N);
    for k = 1 : N
        Pk = P(:, :, k+1);
        K(:, :, k) = -(R+B'*Pk*B)\(B'*Pk*A);
    end

    % compute state x_k and control u_k forward
    u = zeros(Nu, N);
    x = zeros(Nx, N+1);
    x(:, 1) = x0;
    for k = 1 : N
        Kk = K(:, :, k);
        u(:, k) = Kk * x(:, k);
        % projection on input constraints [umin,umax]
        if nargin == 9
            u(:, k) = max(min(u(:,k),umax),umin);
        end
        x(:, k+1) = A * x(:, k) + B * u(:, k);
    end
end
