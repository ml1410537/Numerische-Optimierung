function out = problem_dfdu_vec(x, u, vec, p)
%PROBLEM_DFDU_VEC Product of transposed Jacobi matrix and vector (df / du)' * vec
% 	out = [0 0 0;
%            1 0 0;
%            0 0 0;
%            0 1 0;
%            0 0 0;
%            0 0 1;
%            0 0 0;
%            -1/(x(5)*cos(x(9)))*cos(x(7)) 0 0;
%            0 0 0;
%            1/x(5)*sin(x(7))*sin(x(9)) 1/x(5)*cos(x(9)) 0];
    out(1) = vec(2) + vec(8) * (-1/(x(5)*cos(x(9)))*cos(x(7))) + vec(10) * (1/x(5)*sin(x(7))*sin(x(9)));
    out(2) = vec(4) + vec(10) * (1/x(5)*cos(x(9)));
    out(3) = vec(6);
    out = [out(1) out(2) out(3)];
end
