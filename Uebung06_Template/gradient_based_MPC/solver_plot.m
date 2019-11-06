function solver_plot(vec, sol, tk)
%SOLVER_PLOT Visualize trajectories
    
    subplot(2, 2, 1);
    plot(vec.t, vec.x); xlabel('t'), ylabel('x');
    hold on
    plot(tk + sol.t, sol.x, 'LineStyle', ':');
    hold off;
    title('States');
    
    subplot(2, 2, 2);
    plot(vec.t, vec.u); xlabel('t'), ylabel('u');
    hold on
    plot(tk + sol.t, sol.u, 'LineStyle', ':');
    hold off;
    title('Controls');
    
    subplot(2, 2, 3);
    plot(vec.t, vec.adj); xlabel('t'), ylabel('adj');
    hold on
    plot(tk + sol.t, sol.adj, 'LineStyle', ':');
    hold off;
    title('Adjoint states');
    
    subplot(2, 2, 4);
    plot(vec.t, vec.J); xlabel('t'), ylabel('J');
    title('Cost function');
end