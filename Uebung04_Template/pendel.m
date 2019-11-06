function vec = pendel(xdelta)
%PENDEL Aufschwingen des inversen Pendels mittels Zwei-Freiheitsgrade-Regelung
%
% Vorlesung 'Numerische Optimierung und modellpr?diktive Regelung'
%
% Template f¨¹r 2-FHG-Regelung
%
% Erg?nzen Sie die '???'-Eintr?ge
%

    if nargin == 0
        xdelta = [0; 0; 0; 0];
    end

    % Pendelparameter
    p.g = 9.81;
    p.l = 1.0;
    p.n = 4;
    p.steps = 200;
    
    % Steuerung / Referenzgr??en
    p.tf = 2;
    vec.t = linspace(0, p.tf, p.steps);
    % Eingangstrajektorie
    c = [-135.128, 625.696, -888.598, 497.538, -96.3861];
    vec.uS = c(1)*vec.t+c(2)*(vec.t).^2+c(3)*(vec.t).^3+c(4)*(vec.t).^4+c(5)*(vec.t).^5;
    % Integration der Systemgleichungen
    p.x0 = [0; 0; -pi; 0];
    [~, xsim] = ode45(@simulation, vec.t, p.x0, [], vec, p, 0);
    vec.xS = xsim';
    
    % Gewichtungsmatrizen
    p.S = diag([1, 1, 1, 1]);
    p.Q = diag([1, 1, 1, 1]);
    p.R = 0.01;
    % Matrix -> Vektor
    Svec = reshape(p.S, p.n^2, 1);
    % Integration der Riccati-DGL in R¨¹ckw?rtszeit
    [~, Psim] = ode45(@riccati_dgl, vec.t(end:-1:1), Svec, [], vec, p);
    vec.P = Psim(end:-1:1, :)';
    
    % Zeitvariante Reglerverst?rkungen
    for i = 1 : p.steps
        [~, B] = linearisierung(vec.t(i), vec, p);
        % Vektor -> Matrix
        P = reshape(vec.P(:, i), p.n, p.n);
        vec.k(:, i) = (p.R)^(-1)*B'*P;
    end

    % Simulation mit Anfangsst?rung
    x0 = p.x0 + xdelta;
    % Offener Kreis
    [~, x_ok] = ode45(@simulation, vec.t, x0, [], vec, p, 0);
    vec.x_ok = x_ok';
    % Geschlossener Kreis
    [~, x_gk] = ode45(@simulation, vec.t, x0, [], vec, p, 1);
    vec.x_gk = x_gk';
    
    % Berechnung des Regelungsanteils
    vec.du = zeros(1, p.steps);
    for i = 1 : p.steps
        vec.du(i) = vec.k(:, i)' * (vec.xS(:, i) - vec.x_gk(:, i));
    end
end

function f = simulation(t, x, vec, p, regelung)
    % Steuerung
    uS = interp1(vec.t, vec.uS, t);
    % Simulation im offenen Kreis
    if regelung == 0
        u = uS;
    % Simulation im geschlossenen Kreis
    else
        xS = interp1(vec.t, vec.xS', t)';
        k = interp1(vec.t, vec.k', t);
        u = uS + (-k*(x-xS));
    end
    % Nichtlineares Pendelmodell
    f = [x(2);u;x(4);3/(2*p.l)*(p.g*sin(x(3))-cos(x(3))*u)];
end

function Pvec = riccati_dgl(t, Pvec, vec, p)
    % Zeitvariante Matrizen A und B
    [A, B] = linearisierung(t, vec, p);
    % Vektor -> Matrix
    P = reshape(Pvec, p.n, p.n);
    % Riccati-Differentialgleichung
    P = P + (- P*A - A'*P + P*B*(p.R)^(-1)*B'*P - p.Q);
    % Matrix -> Vektor
    Pvec = reshape(P, p.n^2, 1);
end

function [A, B] = linearisierung(t, vec, p)
    % Momentane Referenzgr??en
    xS = interp1(vec.t, vec.xS', t)';
    uS = interp1(vec.t, vec.uS, t);
    % Zeitvariante Matrizen A und B
    A = [0 1 0 0;0 0 0 0;0 0 0 1;0 0 3/(2*p.l)*(p.g*cos(xS(3))+sin(xS(3))*uS) 0];
    B = [0;1;0;-3/(2*p.l)*cos(xS(3))];
end
