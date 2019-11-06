function vec = schiff_read_AMPL(file)
% ----------------------------------------------------------------
    if ~exist('file', 'var') || isempty(file)
        file = 'schiff_AMPL_sol.txt';
    end

    % Schiffsparameter
    p = schiff_param;

    % Daten lesen und in Struktur vec speichern
    sol = load(file);
    vec.t = sol(:, 1)';
    vec.x = sol(:, [2,3,4])';
    vec.u = sol(:, end)';

    % Trajektorien darstellen
    schiff_plot(vec, p)
end