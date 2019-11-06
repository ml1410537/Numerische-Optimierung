function plot_pendel(vec)
%PLOT_PENDEL Visualisierung des inversen Pendels

    % Trajektorien
    cc = get(gca, 'colororder');
    figure(1);
    
    subplot(2, 3, 1);
    plot(vec.t, vec.xS(1,:), 'k-');
    hold on
    plot(vec.t, vec.x_ok(1,:), '--', 'Color', cc(1,:));
    plot(vec.t, vec.x_gk(1,:), '--', 'Color', cc(2,:));
    hold off
    xlabel('Zeit t'); title('Zustand y');
    
    subplot(2, 3, 2);
    plot(vec.t, vec.xS(2,:), 'k-');
    hold on
    plot(vec.t, vec.x_ok(2,:), '--', 'Color', cc(1,:));
    plot(vec.t, vec.x_gk(2,:), '--', 'Color', cc(2,:));
    hold off
    xlabel('Zeit t'); title('Zustand \dot y');
    
    subplot(2, 3, 3);
    plot(vec.t, vec.uS, 'k-');
    hold on
    plot(vec.t, vec.uS, '--', 'Color', cc(1,:));
    plot(vec.t, vec.uS + vec.du, '--', 'Color', cc(2,:));
    hold off
    xlabel('Zeit t'); title('Stellgröße u');
    
    subplot(2, 3, 4);
    plot(vec.t, vec.xS(3,:), 'k-');
    hold on
    plot(vec.t, vec.x_ok(3,:), '--', 'Color', cc(1,:));
    plot(vec.t, vec.x_gk(3,:), '--', 'Color', cc(2,:));
    hold off
    xlabel('Zeit t'); title('Zustand \phi');
    
    subplot(2, 3, 5);
    plot(vec.t, vec.xS(4,:), 'k-');
    hold on
    plot(vec.t, vec.x_ok(4,:), '--', 'Color', cc(1,:));
    plot(vec.t, vec.x_gk(4,:), '--', 'Color', cc(2,:));
    hold off
    xlabel('Zeit t'); title('Zustand \dot \phi');

    subplot(2, 3, 6);
    plot(vec.t, vec.k);
    xlabel('Zeit t'); title('Verstärkung k');
    
    % Pendelbewegung
    for i = 1 : length(vec.t)
        figure(2);
        clf;
        plot([-2, 2], [-0.05, -0.05], 'k--'); 
        hold on;
        plot_step(vec.xS(:, i), 'k');
        plot_step(vec.x_ok(:, i), cc(1,:));
        plot_step(vec.x_gk(:, i), cc(2,:));
        axis equal
        xlim([-2 2]);
        ylim([-1.1 1.1]);
    end
end

function plot_step(x, color)
    plot([x(1) x(1)+sin(x(3))], [0 cos(x(3))], '-', 'Color', color);
    rectangle('Position', [x(1)-0.1 -0.05 0.2 0.1], 'EdgeColor', color);
end

