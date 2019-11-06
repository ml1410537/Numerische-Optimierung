function out = problem_Vfct(x, p)
%PROBLEM_VFCT Terminal cost V(x(T))
	out = 1/2*(x-p.xdes)'*p.S*(x-p.xdes);
end
