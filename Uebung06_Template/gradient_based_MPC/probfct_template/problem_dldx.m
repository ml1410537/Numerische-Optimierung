function out = problem_dldx(x, u, t, p)
%PROBLEM_DLDX Partial derivatives of integral cost dl / dx
	out = p.Q*(x-p.xdes);
end
