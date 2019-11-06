function out = solver_gradient(tvec, xvec, uvec, adjvec, p)
%SOLVER_GRADIENT Compute gradient of Hamilton function
    Nhor = length(tvec);
    out = zeros(p.Nu, Nhor);
    for k = 1 : Nhor
        out(:,k) = problem_dldu(xvec(:,k), uvec(:,k), tvec(k), p) + ...
                   problem_dfdu_vec(xvec(:,k), uvec(:,k), adjvec(:,k), p);
    end
end