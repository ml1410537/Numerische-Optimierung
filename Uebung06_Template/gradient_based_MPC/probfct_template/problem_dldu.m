function out = problem_dldu(x, u, t, p)
%PROBLEM_DLDU Partial derivatives of integral cost dl / du
	out = p.R*u;
end
