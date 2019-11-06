%
% Vorlesung 'Numerische Optimierung und modellprädiktive Regelung'
%
% Template zu Aufgabe 2.4 und 2.5
%
% Ergänzen Sie die '???'-Einträge
%
clear variables;
close all;
clc;

% input
names = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L'};
activities = 1 : length(names);
durations = [1, 3, 2, 4, 4, 2, 3, 4, 1, 2, 3, 2];
predecessors = [1, 2; % A,B
                1, 3; % A,C
                1, 4; % A,D
                2, 5; % B,E
                3, 5; % C,E
                3, 6; % C,F
                4, 6; % D,F
                5, 7; % E,G
                5, 8; % E,H
                6, 8; % F,H
                6, 9; % F,I
                7, 10; % G,J
                8, 10; % H,J
                8, 11; % H,K
                9, 11; % I,K
                10, 12; % J,L
                11, 12]; % K,L
            
% validate input
for i = 1 : length(names)
    fprintf('Activity: %s\t Duration: %i\t Predecessors: ', names{i}, durations(i));
    for j = 1 : size(predecessors, 1)
        if predecessors(j, 2) == i
            fprintf('%s,', names{predecessors(j, 1)});
        end
    end
    fprintf('\n');
end
fprintf('\n');
    
% create linear programming problem

[A, b, c] = cpm_problem(activities, durations, predecessors);

% solve problem
%sol = simplex(A,b,c);
%x0 = ones(length(c),1);
%[A, b, c, x0] = random_LP(100, 100);
%opt.MaxIter = 200;
%opt.BarrierParamDecreaseFactor = 0.001;
%sol = interiorpoint(A,b,c,x0,opt);
sol = simplex(A,b,c);
fprintf('Optimization finished in %i iterations\n\n', sol.iter);


% output
for i = 1 : length(names)
    fprintf('Activity: %s\t Start time: %.2f\t End time: %.2f\n', names{i}, sol.x(i), sol.x(i)+durations(i));
end
fprintf('Total duration: %.2f\n', sol.x(length(names)+1));

