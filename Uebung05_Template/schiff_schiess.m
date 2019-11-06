function [vec, sol, p] = schiff_schiess(umin, umax, adj0)
% ----------------------------------------------------------------
    % Schiffsparameter
    p = schiff_param;
    % Speichere Stellgrö?enbeschränkungen in Parametern (?bergabe in Grad)
    p.umin = umin * pi/180;
    p.umax = umax * pi/180;

    % Optionen für fsolve
    optoptim = optimset('Display', 'iter');
    % Optionen für Integration mit hoher Genauigkeit
    optode   = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
    % Aufruf von fsolve mit Startlösung adj0
    sol      =  fsolve(@schiff_res,adj0,optoptim,p,optode);

    % Speichere optimale Lösung in vec
    vec.t   = linspace(0, p.tf, 100);
    % Letzte Integration für optimale Lösung
    [~, X]  = ode113(@schiff_kangln,vec.t,[p.x0;sol],optode,p);
    vec.x   = X(:,1:3)';
    vec.adj = X(:,4:6)';
    % Berechnung der optimalen Steuerung
    for i = 1 : length(vec.t)
      vec.u(:,i) = schiff_uopt(vec.x(:,i), vec.adj(:,i), p);
    end
end

function res = schiff_res(adj0, p, optode)
% ----------------------------------------------------------------
    % Randbedingungen für fsolve
    X0    = [p.x0; adj0];
    % Integration der kanonischen Gleichungen
    [~, X] =  ode113(@schiff_kangln,[0,p.tf],X0,optode,p);
    xf    = X(end, 1:3)';
    adjf  = X(end, 4:6)';
    % Endbedingungen in Residuenform
    res   = adjf - p.S*(xf-p.xf);
end