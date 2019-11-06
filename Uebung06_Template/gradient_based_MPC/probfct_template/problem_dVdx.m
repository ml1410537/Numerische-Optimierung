function out = problem_dVdx(x, p)
%PROBLEM_DVDX Partial derivatives of terminal cost dV / dx
	out = p.S*(x-p.xdes);
end
