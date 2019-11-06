
a = 10;
X0 = [0.1;0.1];
opt = optimset('Display','iter');
Xopt = fsolve(@fcn,X0,opt,a)

function F = fcn(x,a)
l = 10;
u0 = x(1);
rl = x(2);
F = [(1+u0^2)^2/(4*u0^3)-a/rl;
     1/4*(3/(4*u0^4)+1/u0^2-7/4+log2(u0))-l/rl];
end