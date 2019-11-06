function schiff_plot(vec, p)
% -------------------------------------------------------------

    % Integration der Differentialgleichungen für Schiffsposition
    [~, way] = ode45(@intpath, vec.t, [0;0], [], vec, p);
    vec.x(4:5, :) = way';

    figure(1);
    subplot(2, 2, 1);
    plot(vec.t, vec.x(1,:) * 180/pi);
    hold on
    xlabel('t');
    ylabel('x_1')

    subplot(2, 2, 2);
    plot(vec.t, vec.x(3,:) * 180/pi);
    hold on
    xlabel('t');
    ylabel('x_3')

    subplot(2, 2, 3);
    plot(vec.t, vec.u * 180/pi);
    hold on
    xlabel('t');
    ylabel('u');

    subplot(2, 2, 4)
    plot(vec.x(4,:), vec.x(5,:));
    xlabel('x_4');
    ylabel('x_5');
    
    % Schiffsposition alle 10 Schritte darstellen
    cc = get(gca, 'colororder');
    len = 2;
    hold on
    for i = 1 : 10 : length(vec.t)
      xs = vec.x(4, i);
      ys = vec.x(5, i);
      psis = vec.x(1, i);
      plot(xs+len*sin(psis)*[-1 1], ys+len*cos(psis)*[-1 1], 'Color', cc(1,:), 'LineWidth',3);
    end
    axis equal
end

function F = intpath(t, x, vec, p)
% -------------------------------------------------------------
    % Zustand zum Zeitpunkt t interpolieren
    xref = interp1(vec.t, vec.x', t)';
    % Schiffsposition (x_4, x_5) berechnen
    F =  [p.v*sin(xref(1)-xref(3));
          p.v*cos(xref(1)-xref(3))];
end
