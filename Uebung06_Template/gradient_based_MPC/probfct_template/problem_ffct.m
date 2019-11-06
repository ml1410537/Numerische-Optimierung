function out = problem_ffct(x, u, p)
%PROBLEM_FFCT Right-hand side of differential equation \dot x = f(x, u)	
	out = [x(2);
           u(1);
           x(4);
           u(2);
           x(6);
           u(3);
           x(8);
           -1/(x(5)*cos(x(9)))*(9.8*sin(x(7))+u(1)*cos(x(7))+2*x(8)*(x(6)*cos(x(9))-x(5)*x(10)*sin(x(9))));
           x(10);
           -1/x(5)*(2*x(6)*x(10)+sin(x(9))*(9.8*cos(x(7))-u(1)*sin(x(7)))-u(1)*sin(x(7))-u(2)*cos(x(9))-x(5)*x(8)^2*sin(x(9))*cos(x(9)))];
end
