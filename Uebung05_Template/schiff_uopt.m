function u = schiff_uopt(x, adj, p)
% ----------------------------------------------------------------
    % Berechnung der optimalen Stellgrö?e ohne Beschränkungen
    u0 =  -1/p.r*[0 p.c(2) 0]*adj;

    % Betrachtung der Beschränkungen
    if u0> p.umin & u0<p.umax
        u = u0;
    elseif u0<=p.umin
        u = p.umin;
    else
        u = p.umax;
    end
end
