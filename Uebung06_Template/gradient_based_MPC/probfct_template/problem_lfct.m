function out = problem_lfct(x, u, t, p)
%PROBLEM_LFCT Integral cost l(x, u, t)
	out = 1/2*((x-p.xdes)'*p.Q*(x-p.xdes)+u'*p.R*u);
end
