function out = solver_project(u, p)
%SOLVER_PROJECT Project control u on feasible set [p.umin, p.umax]
    out = u;
    Nhor = size(u, 2);
    for k = 1 : Nhor
        out(:,k) = max(p.umin, min(u(:,k), p.umax));
    end
end