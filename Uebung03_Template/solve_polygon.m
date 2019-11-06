function [xsol, fsol, mult] = solve_polygon(N)
%SOLVE_POLYGON Solve largest small polygon problem with N edges

% compute initial polygon
angle0 = linspace(2*pi/N, 2*pi, N);
xx0 = 0.5 * sin(angle0);
yy0 = 0.5 - 0.5 * cos(angle0);

r0 = sqrt(xx0.^2 + yy0.^2);
th0 = zeros(1, N);

for i = 1 : N
  if xx0(i) > 0
    th0(i) = atan(yy0(i) / xx0(i));
  elseif xx0(i) == 0
    th0(i) = pi/2;
  else
    th0(i) = pi - atan(-yy0(i) / xx0(i));
  end
end


% plot
figure(1); 
clf;
plot_polygon(r0, th0, 'r--');
hold on;
rectangle('Position', [-0.5 0 1 1], 'Curvature', [1 1], 'EdgeColor', 'k');

% create nonlinear constraints
[C, ~] = nonlcon([r0, th0], N, 0);
Nc = length(C);

% solve using fmincon
opt = optimset('Display', 'iter', 'TolX', 1e-12, 'TolFun', 1e-12, 'TolCon', 1e-12, 'Algorithm', 'active-set');
[xsol, fsol, ~, ~, mult] = ...
    fmincon(@cost, [r0, th0], [], [], [], [], [], [], @nonlcon, opt, N, Nc);

% plot
r  = xsol(1:N);
th = xsol(N+1:end);
plot_polygon(r, th, 'b');
title(['Largest small polygon, N = ', int2str(length(r))]);

fprintf('Optimaler Kostenwert: %g\n', fsol);

% ---------------------------------------------------
function plot_polygon(r, th, mycolor)
% coordinates
xx = r .* cos(th);
xx(end+1) = xx(1);
yy = r .* sin(th); 
yy(end+1) = yy(1);
plot(xx, yy, mycolor);
hold on;
axis equal;

% ---------------------------------------------------
function f = cost(x, N, Nc)
r = x(1:N);
th = x(N+1:end);
f = 0;
for i=1:N-1
   f = f - 0.5*r(i)*r(i+1)*sin(th(i+1)-th(i));
end

% ---------------------------------------------------
function [C, Ceq] = nonlcon(x, N, Nc)
r = x(1:N);
th = x(N+1:end);

% equality constraints
Ceq = [r(end);th(end)-pi];

if Nc > 0
    C = zeros(Nc,1); 
end

% inequality constraints
k = 1;
for i = 1 : N-2
  for j = i+1 : N-1
    C(k) = r(i)^2 + r(j)^2 - 2*r(i)*r(j)*cos(th(j)-th(i))-1;
    k = k+1;
  end
end

for i = 1 : N-2
  C(k) = th(i)-th(i+1);
  k = k+1;
end

for i = 1 : N-1
  C(k) = -th(i);
  C(k+1) = th(i)-pi;
  k = k+2;
end

for i = 1 : N-1
  C(k) = -r(i);
  C(k+1) = r(i)-1;
  k = k+2;
end
