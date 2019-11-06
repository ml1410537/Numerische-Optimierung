function f = schiff_kangln(t, X, p)
% ----------------------------------------------------------------
    x = X(1:3);
    adj = X(4:6);
    % Berechnung der optimalen Steuerung
    u = schiff_uopt(x, adj, p);
    % Kanonische Gleichungen
    f =  [ fsys(x,u,p); 
           -p.Q*(x-p.xf)-dfdx(x,adj,u,p)'*adj ];
end

function f = fsys(x, u, p)
% ----------------------------------------------------------------
    % Berechnung der Schiffsdynamik
    f = [ x(2); 
          p.c(1)*x(2)+p.c(2)*u;
          p.c(3)*p.v*abs(x(3))*x(3)+p.c(4)*x(2)];
end

function J = dfdx(x, adj, u, p)
% ----------------------------------------------------------------
    % Berechnung der Jacobi-Matrix der Schiffsdynamik
    J = [ 0, 1, 0; 
          0, p.c(1), 0;
          0, p.c(4), 2*p.c(3)*p.v*x(3)*sign(x(3)) ];
end
