function [vec, sol, p] = schiff_rwa(umin, umax)
% ----------------------------------------------------------------
    % Schiffsparameter
    p = schiff_param;
    % Speichere Stellgrö?enbeschränkungen in Parametern (?bergabe in Grad)
    p.umin = umin * pi/180;
    p.umax = umax * pi/180;

    % Optionen für bvp4c
    opt  = bvpset('Stats', 'on', 'RelTol', 1e-4);
    % Startlösung für Zustände und adjungierte Zustände
    sol0 = bvpinit(linspace(0, p.tf, 30), [p.x0; zeros(3,1)]);
    % Aufruf von bvp4c
    sol  =  bvp4c(@schiff_kangln, @schiff_rb, sol0, opt, p);

    % Speichere optimale Lösung in vec
    vec.t   = linspace(0, p.tf, 100);
    vec.x   = interp1(sol.x, sol.y(1:3,:)', vec.t)';
    vec.adj = interp1(sol.x, sol.y(4:6,:)', vec.t)';
    % Berechnung der optimalen Steuerung
    for i = 1 : length(vec.t)
      vec.u(:,i) = schiff_uopt(vec.x(:,i), vec.adj(:,i), p);
    end
end

function res = schiff_rb(X0, Xf, p)
% ----------------------------------------------------------------
    % Randbedingungen für bvp4c
    x0 = X0(1:3);
    xf = Xf(1:3);
    adjf = Xf(4:6);
    % Residuum
    res =  [ x0 - p.x0;
            adjf - p.S*(xf-p.xf)];
end
