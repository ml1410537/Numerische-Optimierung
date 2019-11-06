function p = schiff_param
% ----------------------------------------------------------------
    % Parameter der Schiffsdynamik
    p.c  = [-0.26, 0.2, -1.87, 0.6];
    % Geschwindigkeit
    p.v  = 3.5;
    % Endzeit (Optimierungshorizont)
    p.tf = 15;

    % Anfangsbedingungen
    p.x0 = [0; 0; 0];
    % (Gewünschte) Endbedingungen
    p.xf = [45*pi/180; 0; 0];

    % Gewichtungen des Endzustands                       
    p.S  = diag([1, 1, 1]);
    % Gewichtungen der Zustände im Integralanteil
    p.Q  = diag([1, 1, 1]);
    % Gewichtung der Stellgröße im Integralanteil
    p.r  = 0.1;
end